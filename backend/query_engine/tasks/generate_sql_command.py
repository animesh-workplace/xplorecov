from celery import shared_task
from vanna.ollama import Ollama
from django.conf import settings
from vanna.chromadb import ChromaDB_VectorStore

class MyVanna(ChromaDB_VectorStore, Ollama):
    def __init__(self, config=None):
        ChromaDB_VectorStore.__init__(self, config=config)
        Ollama.__init__(self, config=config)

@shared_task
def query_using_ai_model(question):
    database_path = str(settings.BASE_DIR / 'database')
    sqlite_database_path = str(settings.BASE_DIR / 'database' / 'db.sqlite3')
    model = MyVanna(config={'model': 'granite3-dense', "path": database_path })
    try:
        print(model.generate_questions())
        sql_query = model.generate_sql(question)
        if(model.is_sql_valid(sql_query)):
            return {"query": sql_query}
        else:
            return {"message": "Query incorrect"}
    except Exception as e:
        return {"exception_occured": e}


