from dotenv import load_dotenv
from celery import shared_task
from vanna.ollama import Ollama
from django.conf import settings
from vanna.chromadb import ChromaDB_VectorStore


class MyVanna(ChromaDB_VectorStore, Ollama):
    def __init__(self, config=None):
        ChromaDB_VectorStore.__init__(self, config=config)
        Ollama.__init__(self, config=config)

@shared_task
def train_ai_model():
    database_path = str(settings.BASE_DIR / 'database')
    sqlite_database_path = str(settings.BASE_DIR / 'database' / 'db.sqlite3')
    model = MyVanna(config={'model': 'granite3-dense', "path": database_path })
    model.connect_to_sqlite(sqlite_database_path)
    df_ddl = model.run_sql("SELECT type, sql FROM sqlite_master WHERE name LIKE 'query_engine%' AND sql is not null")
    for ddl in df_ddl['sql'].to_list():
        model.train(ddl=ddl)

    # VirusStrain Table
    model.train(ddl="""
        CREATE TABLE IF NOT EXISTS query_engine_virusstrain (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strain VARCHAR(255) NOT NULL,
            gisaid_epi_isl VARCHAR(50) UNIQUE NOT NULL,
            genbank_accession VARCHAR(50),
            collection_date DATE,
            host VARCHAR(100),
            age INTEGER,
            sex VARCHAR(10),
            coverage DECIMAL(10, 5) CHECK (coverage BETWEEN 0 AND 1),
            nextclade_clade VARCHAR(50),
            nextclade_lineage VARCHAR(50),
            pangolin_lineage VARCHAR(50),
            pangolin_qc_status VARCHAR(20) CHECK (qc_status IN ('fail', 'pass'))
        );
    """)

    # VirusStrain Indexes
    model.train(ddl="CREATE INDEX IF NOT EXISTS virusstrain_collection_date_idx ON query_engine_virusstrain (collection_date);")
    model.train(ddl="CREATE INDEX IF NOT EXISTS virusstrain_nextclade_clade_idx ON query_engine_virusstrain (nextclade_clade);")
    model.train(ddl="CREATE INDEX IF NOT EXISTS virusstrain_nextclade_lineage_idx ON query_engine_virusstrain (nextclade_lineage);")
    model.train(ddl="CREATE INDEX IF NOT EXISTS virusstrain_pangolin_lineage_idx ON query_engine_virusstrain (pangolin_lineage);")

    # Location Table
    model.train(ddl="""
        CREATE TABLE IF NOT EXISTS query_engine_location (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strain_id INTEGER UNIQUE REFERENCES query_engine_virusstrain(id) ON DELETE CASCADE,
            continent VARCHAR(100),
            country VARCHAR(100),
            state VARCHAR(100),
            district VARCHAR(255),
            continent_exposure VARCHAR(100),
            country_exposure VARCHAR(100),
            state_exposure VARCHAR(100)
        );
    """)

    # Location Indexes
    model.train(ddl="CREATE INDEX IF NOT EXISTS location_state_idx ON query_engine_location (state);")
    model.train(ddl="CREATE INDEX IF NOT EXISTS location_country_idx ON query_engine_location (country);")
    model.train(ddl="CREATE INDEX IF NOT EXISTS location_district_idx ON query_engine_location (district);")

    # AAMutation Table
    model.train(ddl="""
        CREATE TABLE IF NOT EXISTS query_engine_aamutation (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strain_id INTEGER REFERENCES query_engine_virusstrain(id) ON DELETE CASCADE,
            protein VARCHAR(50) CHECK (protein IN (
                'N', 'S', 'E', 'M', 'ORF1a', 'ORF1b', 'ORF3a', 'ORF3b', 'ORF6', 'ORF7a', 
                'ORF7b', 'ORF8', 'ORF9b', 'ORF10', 'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 
                'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP12', 'NSP13', 'NSP14', 
                'NSP15', 'NSP16'
            )),
            reference_aa VARCHAR(10),
            position INTEGER,
            altered_aa VARCHAR(10)
        );
    """)

    # AAMutation Index
    model.train(ddl="CREATE INDEX IF NOT EXISTS aamutation_protein_position_ref_alt_idx ON query_engine_aamutation (protein, position, reference_aa, altered_aa);")

    # NucMutation Table
    model.train(ddl="""
        CREATE TABLE IF NOT EXISTS query_engine_nucmutation (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strain_id INTEGER REFERENCES query_engine_virusstrain(id) ON DELETE CASCADE,
            reference_base CHAR(1),
            position INTEGER,
            altered_base CHAR(1)
        );
    """)

    # NucMutation Index
    model.train(ddl="CREATE INDEX IF NOT EXISTS nucmutation_ref_pos_alt_idx ON query_engine_nucmutation (reference_base, position, altered_base);")

    # StructuralVariation Table
    model.train(ddl="""
        CREATE TABLE IF NOT EXISTS query_engine_structuralvariation (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strain_id INTEGER UNIQUE REFERENCES query_engine_virusstrain(id) ON DELETE CASCADE,
            total_insertions INTEGER DEFAULT 0,
            insertion_positions TEXT,
            total_frameshifts INTEGER DEFAULT 0,
            frameshift_positions TEXT,
            total_missing INTEGER DEFAULT 0,
            missing_ranges TEXT
        );
    """)

    # StructuralVariation Index
    model.train(ddl="CREATE INDEX IF NOT EXISTS structuralvariation_insertions_frameshifts_missing_idx ON query_engine_structuralvariation (total_insertions, total_frameshifts, total_missing);")

    # QualityMetrics Table
    model.train(ddl="""
        CREATE TABLE IF NOT EXISTS query_engine_qualitymetrics (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strain_id INTEGER UNIQUE REFERENCES query_engine_virusstrain(id) ON DELETE CASCADE,
            total_non_acgtns INTEGER DEFAULT 0,
            non_acgtn_positions TEXT,
            partially_aliased BOOLEAN DEFAULT 0,
            scorpio_call VARCHAR(50),
            scorpio_support DECIMAL(5, 2),
            scorpio_conflict DECIMAL(5, 2)
        );
    """)

    # Documentation on key business terminology
    # VirusStrain Table Documentation
    model.train(documentation="""
    The 'query_engine_virusstrain' table is designed to store essential metadata for individual viral strains. 
    Each strain has a unique identifier, 'gisaid_epi_isl', which ensures its distinct identity in the database. 
    Other columns provide detailed information, such as:
    - 'strain': The name of the strain, limited to 255 characters.
    - 'gisaid_epi_isl': A unique, mandatory identifier from the GISAID database, with a maximum length of 50 characters.
    - 'genbank_accession': An optional field for the GenBank accession number.
    - 'collection_date': The date when the sample was collected which contains the date in the YYYY-MM-DD format.
    - 'host': The host organism from which the sample was collected, with a limit of 100 characters.
    - 'age': The age of the host, if available.
    - 'sex': The sex of the host, with a 10-character limit.
    - 'coverage': A decimal value representing the sequencing coverage, constrained between 0 and 1, with up to 5 decimal places.
    - 'nextclade_clade': A fields for the Nextclade clade assignment, each with a maximum length of 50 characters.
    - 'nextclade_lineage': A fields for the Nextclade lineage assignments, each with a maximum length of 50 characters.
    - 'pangolin_lineage': An field for Pangolin lineage classification, up to 50 characters.
    - 'pangolin_qc_status': A quality control performed by pangolin assigning tool status field, with choices of 'pass' or 'fail'.
    - 'nextclade_qc_status': A quality control performed by nextclade clade/lineage assigning tool status with choices of 'good', 'mediocre', or 'bad'.
    Indexes are created on 'collection_date', 'nextclade_clade', 'nextclade_lineage', and 'pangolin_lineage' 
    to enable efficient filtering based on these fields.
    """)

    # Location Table Documentation
    model.train(documentation="""
    The 'query_engine_location' table stores detailed geographic metadata linked to each viral strain, with a 
    one-to-one relationship to the 'query_engine_virusstrain' table. This table includes:
    - 'strain_id': A unique reference to a virus strain in the 'query_engine_virusstrain' table.
    - 'continent': The continent where the strain originated, up to 100 characters.
    - 'country': The country of origin, limited to 100 characters.
    - 'state': The state or region of origin, limited to 100 characters.
    - 'district': The specific district where the sample was collected, up to 255 characters.
    - 'continent_exposure', 'country_exposure', 'state_exposure': Optional fields that represent the 
    continent, country, and state where exposure is believed to have occurred.
    Indexes are created on 'state', 'country', and 'district' for faster querying and filtering based on 
    geographic fields.
    """)

    # AAMutation Table Documentation
    model.train(documentation="""
    The 'query_engine_aamutation' table records amino acid mutations linked to each viral strain. It includes 
    information on the specific protein, position, and amino acid change. Key columns are:
    - 'strain_id': A reference to the related strain in 'query_engine_virusstrain'.
    - 'protein': The protein where the mutation occurs, with predefined choices including:
        - 'N' for Nucleocapsid protein,
        - 'S' for Spike protein,
        - 'E' for Envelope protein,
        - 'M' for Membrane protein,
        - 'ORF1a' for ORF1a Polyprotein,
        - 'ORF1b' for ORF1b Polyprotein,
        - 'ORF3a' for ORF3a Protein,
        - 'ORF3b' for ORF3b Protein,
        - 'ORF6' for ORF6 Protein,
        - 'ORF7a' for ORF7a Protein,
        - 'ORF7b' for ORF7b Protein,
        - 'ORF8' for ORF8 Protein,
        - 'ORF9b' for ORF9b Protein,
        - 'ORF10' for ORF10 Protein,
        - 'NSP1' for Non-structural Protein 1,
        - 'NSP2' for Non-structural Protein 2,
        - 'NSP3' for Non-structural Protein 3,
        - 'NSP4' for Non-structural Protein 4,
        - 'NSP5' for Main Protease (NSP5),
        - 'NSP6' for Non-structural Protein 6,
        - 'NSP7' for Non-structural Protein 7,
        - 'NSP8' for Non-structural Protein 8,
        - 'NSP9' for Non-structural Protein 9,
        - 'NSP10' for Non-structural Protein 10,
        - 'NSP12' for RNA-dependent RNA Polymerase (NSP12),
        - 'NSP13' for Helicase (NSP13),
        - 'NSP14' for 3'-to-5' Exonuclease (NSP14),
        - 'NSP15' for Endoribonuclease (NSP15), and
        - 'NSP16' for 2'-O-Methyltransferase (NSP16).
    - 'reference_aa': The original amino acid before mutation, limited to 10 characters.
    - 'position': The position in the protein sequence where the mutation occurs.
    - 'altered_aa': The mutated amino acid, up to 10 characters.
    An index on 'protein', 'position', 'reference_aa', and 'altered_aa' is included to facilitate efficient 
    searches for specific mutations.
    """)


    # NucMutation Table Documentation
    model.train(documentation="""
    The 'query_engine_nucmutation' table stores nucleotide mutations associated with each viral strain. Each 
    mutation entry consists of:
    - 'strain_id': A reference to the viral strain in the 'query_engine_virusstrain' table.
    - 'reference_base': The original nucleotide base before mutation, a single character.
    - 'position': The genomic position of the mutation.
    - 'altered_base': The nucleotide base after mutation, a single character.
    An index on 'reference_base', 'position', and 'altered_base' enables fast lookups for nucleotide mutation 
    queries.
    """)

    # StructuralVariation Table Documentation
    model.train(documentation="""
    The 'query_engine_structuralvariation' table records structural variation data for each viral strain, 
    including counts of insertions, frameshifts, and missing segments. This table contains:
    - 'strain_id': A unique reference to the viral strain in 'query_engine_virusstrain'.
    - 'total_insertions': An integer count of total insertions found in the genome.
    - 'insertion_positions': A JSON-formatted string representing insertion positions.
    - 'total_frameshifts': An integer count of total frameshifts.
    - 'frameshift_positions': A JSON-formatted string representing frameshift positions.
    - 'total_missing': An integer count of missing segments.
    - 'missing_ranges': A JSON-formatted string listing the ranges of missing segments.
    Indexes on 'total_insertions', 'total_frameshifts', and 'total_missing' allow for efficient querying of 
    structural variation statistics.
    """)

    # QualityMetrics Table Documentation
    model.train(documentation="""
    The 'query_engine_qualitymetrics' table provides quality metrics for each strain, capturing sequence quality 
    and additional evaluation data. This includes:
    - 'strain_id': A unique reference to the virus strain in the 'query_engine_virusstrain' table.
    - 'total_non_acgtns': An integer representing the count of non-ACGTN characters in the sequence.
    - 'non_acgtn_positions': A JSON-formatted string listing positions of non-ACGTN characters.
    - 'partially_aliased': A boolean indicating whether the sequence is partially aliased.
    - 'scorpio_call': Classification by the Scorpio tool, limited to 50 characters.
    - 'scorpio_support' and 'scorpio_conflict': Decimal fields for the support and conflict scores 
    associated with the Scorpio call.
    This table provides detailed quality assessment information but does not include indexes as it serves 
    primarily for quality tracking purposes.
    """)

    model.train(documentation="""When lineage is mentioned use the columns nextclade_lineage and pangolin_lineage from query_engine_virusstrain""")

    # SQl based training on questions
    model.train(question="get all the sequences with XBB lineage between the year 2023 and 2024 of the state West Bengal", sql="SELECT * FROM query_engine_virusstrain AS vs JOIN query_engine_location AS loc ON vs.id = loc.strain_id WHERE vs.nextclade_lineage = 'XBB' AND vs.pangolin_lineage = 'XBB' AND vs.collection_date BETWEEN '2023-01-01' AND '2024-12-31' AND loc.country = 'India' AND loc.state = 'West Bengal'")
    # model.train(question="find the sequences with that qc failed but the lineage was assigned", sql="SELECT * FROM query_engine_virusstrain AS vs WHERE vs.qc_status = 'bad' AND vs.lineage IS NOT NULL")

    return { "message": "Training complete" }
