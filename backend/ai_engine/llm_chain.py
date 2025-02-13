import json
from openai import OpenAI
from guardrails import Guard
from dataclasses import dataclass
from datetime import datetime, date
from django.db.models import Model, QuerySet
from guardrails.errors import ValidationError
from django.core.serializers import serialize
from typing import Dict, List, Optional, Union
from guardrails.hub import NSFWText, WebSanitization, DetectJailbreak
from query_engine.models import *
from django.db.models import (
    F,
    Q,
    Avg,
    Max,
    Min,
    Sum,
    Count,
    StdDev,
    Variance,
    Subquery,
    OuterRef,
    Exists,
)


@dataclass
class AgentResponse:
    success: bool
    message: str
    data: Optional[Dict] = None


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (datetime, date)):
            return obj.isoformat()
        return super().default(obj)


class GuardAgent:
    def validate_request(self, user_request: str) -> AgentResponse:
        """
        Validates if the request is safe and feasible to process.
        Returns TRUE/FALSE with reasoning.
        """
        guard = Guard().use_many(
            DetectJailbreak(),
            WebSanitization(),
            NSFWText(threshold=0.8, validation_method="sentence"),
        )

        try:
            guard.validate(user_request)
            return AgentResponse(success=True, message="Request approved")
        except ValidationError as e:
            return AgentResponse(
                success=False,
                message="""
                It looks like your request contains some content that has been flagged as inappropriate. 
                We strive to keep the content safe and respectful for all users. 
                Could you please rephrase or modify your request?
                """,
            )


class DBPromptEngineer:
    def generate_db_prompt(self, user_request: str, db_schema: Dict) -> str:
        """
        Generates database-specific prompt based on the request and schema.
        """
        system_prompt = f"""
        You are an expert prompt engineering AI assistant specialized in database development. 
        Your role is to enhance user prompts by incorporating database schema details to create comprehensive, technically precise prompts for code generation 
        in Django ORM.
        
        When presented with a user prompt and database schema, you should:
        - Analyze the user's intent and required functionality carefully
        - Identify all relevant tables and columns from the provided schema
        - Consider relationships between tables, data types, and constraints

        Keep the prompt short but descriptive. JUST PROVIDE THE PROMPT AND NOTHING ELSE

        Example Prompt:

        User Prompt - Identify strains that belong to the AY.24 lineage
        Your Prompt - Write a Django ORM query to retrieve all virus strains that belong to the AY.24 lineage, using the 'virusstrain' table and the 'nextclade_lineage' field. Include related 'location' and 'quality_metrics' tables in the query to provide additional strain information. 

        Query the following fields:
        - virusstrain.strain
        - virusstrain.gisaid_epi_isl
        - location.country
        - location.continent
        - quality_metrics.scorpio_support
        - quality_metrics.scorpio_conflict

        Filter the results to include only strains with the 'nextclade_lineage' set to 'AY.24'.

        """

        client = OpenAI(
            base_url="http://10.10.6.80/xplorecov/ai/code/v1",
            api_key="sk-no-key-required",
        )
        completion = client.chat.completions.create(
            temperature=1,
            model="LLaMA_CPP",
            messages=[
                {
                    "role": "system",
                    "content": system_prompt,
                },
                {
                    "role": "user",
                    "content": f"""
                    Given the database schema:
                    {db_schema}
                    Generate the advanced prompt based on the user's request {user_request}
                """,
                },
            ],
        )

        response = completion.choices[0].message.content
        return response


class CoderAgent:
    def generate_orm_code(self, db_prompt: str, db_schema: Dict) -> str:
        """
        Generates Django ORM code based on the engineered prompt.
        """
        client = OpenAI(
            base_url="http://10.10.6.80/xplorecov/ai/code/v1",
            api_key="sk-no-key-required",
        )

        system_prompt = f"""
            You are an expert in Django ORM and database queries. Your task is to generate efficient and optimized Django ORM code for the user. 
            Given the database schema:
            {db_schema} 

            The database in the project is SQLITE. 
            Store the result always in a variable called orm_result.
            JUST PROVIDE THE DJANGO ORM CODE WITHOUT EXPLANATION without code formatting.
        """

        completion = client.chat.completions.create(
            temperature=1,
            model="LLaMA_CPP",
            messages=[
                {
                    "role": "system",
                    "content": system_prompt,
                },
                {
                    "role": "user",
                    "content": db_prompt,
                },
            ],
        )

        response = completion.choices[0].message.content
        return response

    def execute_query(self, orm_code: str) -> QuerySet:
        """
        Executes the generated ORM code safely and returns the results.
        """
        local_vars = {}
        try:
            # Execute the generated ORM code
            exec(orm_code, globals(), local_vars)
            result = local_vars.get("orm_result")
            return result
        except Exception as e:
            raise Exception(f"Query execution failed: {str(e)}")


class SummarizerAgent:
    def summarize_data(self, data: List) -> str:
        """
        Summarizes the retrieved data into a user-friendly format.
        """
        client = OpenAI(
            base_url="http://10.10.6.80/xplorecov/ai/code/v1",
            api_key="sk-no-key-required",
        )

        system_prompt = f"""
            You are a expert data engineer of SARS-COV-2. You task is to summarize the data given to you and provide certain inference on it. 
            Summarize the following dataset for the user and write the output in short paragraph.
        """

        completion = client.chat.completions.create(
            temperature=1,
            model="LLaMA_CPP",
            messages=[
                {"role": "system", "content": system_prompt},
                {
                    "role": "user",
                    "content": json.dumps(
                        list(data), cls=DateTimeEncoder, default=lambda o: float(o)
                    ),
                },
            ],
        )

        response = completion.choices[0].message.content
        return response


class AgentChain:
    def __init__(self):
        self.guard = GuardAgent()
        self.coder = CoderAgent()
        self.summarizer = SummarizerAgent()
        self.prompt_engineer = DBPromptEngineer()

    def process_request(self, user_request: str, db_schema: Dict) -> AgentResponse:
        """
        Orchestrates the entire agent chain process.
        """
        # Step 1: Validate request
        validation = self.guard.validate_request(user_request)
        if not validation.success:
            return validation

        try:
            # Step 1: Generate DB-specific prompt
            db_prompt = self.prompt_engineer.generate_db_prompt(user_request, db_schema)
            print(db_prompt)

            # Step 2: Generate ORM code
            orm_code = self.coder.generate_orm_code(db_prompt, db_schema)
            print(orm_code)

            # Step 3: Execute ORM code
            query_results = self.coder.execute_query(orm_code)
            print(query_results)

            # Step 4: Summarize results
            summary = self.summarizer.summarize_data(query_results)

            return AgentResponse(
                success=True,
                message="Request processed successfully",
                data={"summary": summary, "raw_data": query_results},
            )

        except Exception as e:
            return AgentResponse(
                success=False, message=f"Error processing request: {str(e)}"
            )
