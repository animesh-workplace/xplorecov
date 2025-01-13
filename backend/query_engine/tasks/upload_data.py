from tqdm import tqdm
import io, pandas, os, shutil
from celery import shared_task
from ..models import (
    VirusStrain,
    Location,
    AAMutation,
    NucMutation,
    StructuralVariation,
    QualityMetrics,
    parse_aamutation_string,
)


@shared_task
def process_csv_file(file_loc):
    # Convert CSV content back to a DataFrame
    df = pandas.read_feather(file_loc)

    # Process CSV file in chunks
    errors = []
    chunk_size = 1000
    total_records = 0

    for i in tqdm(range(0, len(df), chunk_size)):
        chunk = df[i : i + chunk_size]
        try:
            total_records += _process_chunk(chunk)
        except Exception as e:
            errors.append(f"Error processing records {i} to {i + chunk_size}: {str(e)}")

    # Remove the temporary folder
    shutil.rmtree(os.path.dirname(file_loc), ignore_errors=True)
    return {"total_records": total_records, "errors": errors if errors else None}


def _process_chunk(chunk):
    location_objects = []
    aa_mutation_objects = []
    nuc_mutation_objects = []
    virus_strain_objects = []
    quality_metrics_objects = []
    structural_variation_objects = []

    records_processed = 0

    for _, row in chunk.iterrows():
        try:
            # Create VirusStrain object without saving
            strain = VirusStrain(
                strain=row["strain"],
                gisaid_epi_isl=row["gisaid_epi_isl"],
                genbank_accession=row.get("genbank_accession"),
                collection_date=(
                    pandas.to_datetime(row["date"]).date()
                    if pandas.notna(row["date"])
                    else None
                ),
                host=row.get("host"),
                age=row.get("age") if pandas.notna(row.get("age")) else None,
                sex=row.get("sex"),
                coverage=(
                    float(row["coverage"])
                    if pandas.notna(row.get("coverage"))
                    else None
                ),
                nextclade_clade=row.get("clade"),
                nextclade_lineage=row.get("Nextclade_pango"),
                pangolin_lineage=row.get("lineage"),
                nextclade_qc_status=row.get("qc.overallStatus"),
                pangolin_qc_status=row.get("qc_status"),
            )
            virus_strain_objects.append(strain)

            # Create Location object linked to strain
            location_objects.append(
                Location(
                    strain=strain,
                    continent=row.get("region"),
                    country=row.get("country"),
                    state=row.get("division"),
                    district=row.get("location"),
                    continent_exposure=row.get("region_exposure"),
                    country_exposure=row.get("country_exposure"),
                    state_exposure=row.get("division_exposure"),
                )
            )

            # Process AA mutations
            if pandas.notna(row.get("aaSubstitutions")):
                aa_mutation_objects.extend(
                    _process_aa_mutations(strain, row["aaSubstitutions"])
                )

            # Process nucleotide mutations
            if pandas.notna(row.get("substitutions")):
                nuc_mutation_objects.extend(
                    _process_nuc_mutations(strain, row["substitutions"])
                )

            # Create StructuralVariation object linked to strain
            structural_variation_objects.append(
                StructuralVariation(
                    strain=strain,
                    total_insertions=(
                        int(row["totalInsertions"])
                        if pandas.notna(row.get("totalInsertions"))
                        else 0
                    ),
                    insertion_positions=row.get("insertions"),
                    total_frameshifts=(
                        int(row["totalFrameShifts"])
                        if pandas.notna(row.get("totalFrameShifts"))
                        else 0
                    ),
                    frameshift_positions=row.get("frameShifts"),
                    total_missing=(
                        int(row["totalMissing"])
                        if pandas.notna(row.get("totalMissing"))
                        else 0
                    ),
                    missing_ranges=row.get("missing"),
                )
            )

            # Create QualityMetrics object linked to strain
            quality_metrics_objects.append(
                QualityMetrics(
                    strain=strain,
                    total_non_acgtns=(
                        int(row["totalNonACGTNs"])
                        if pandas.notna(row.get("totalNonACGTNs"))
                        else 0
                    ),
                    non_acgtn_positions=row.get("nonACGTNs"),
                    partially_aliased=(
                        bool(row["partiallyAliased"])
                        if pandas.notna(row.get("partiallyAliased"))
                        else False
                    ),
                    scorpio_call=row.get("scorpio_call"),
                    scorpio_support=(
                        float(row["scorpio_support"])
                        if pandas.notna(row.get("scorpio_support"))
                        else None
                    ),
                    scorpio_conflict=(
                        float(row["scorpio_conflict"])
                        if pandas.notna(row.get("scorpio_conflict"))
                        else None
                    ),
                )
            )

            records_processed += 1

        except Exception as e:
            raise Exception(f"Error processing row {records_processed}: {str(e)}")

    # Bulk insert all objects
    VirusStrain.objects.bulk_create(
        virus_strain_objects, batch_size=len(virus_strain_objects)
    )
    Location.objects.bulk_create(location_objects, batch_size=len(location_objects))
    AAMutation.objects.bulk_create(
        aa_mutation_objects, batch_size=len(aa_mutation_objects)
    )
    NucMutation.objects.bulk_create(
        nuc_mutation_objects, batch_size=len(nuc_mutation_objects)
    )
    StructuralVariation.objects.bulk_create(
        structural_variation_objects, batch_size=len(structural_variation_objects)
    )
    QualityMetrics.objects.bulk_create(
        quality_metrics_objects, batch_size=len(quality_metrics_objects)
    )

    return records_processed


def _process_aa_mutations(strain, mutations_str):
    aa_mutation_objects = []
    mutations = mutations_str.split(",")
    for mutation in mutations:
        mutation_string = parse_aamutation_string(mutation)
        aa_mutation_objects.append(
            AAMutation(
                strain=strain,
                protein=mutation_string["protein"],
                reference_aa=mutation_string["reference"],
                position=mutation_string["position"],
                altered_aa=mutation_string["alternate"],
            )
        )
    return aa_mutation_objects


def _process_nuc_mutations(strain, mutations_str):
    nuc_mutation_objects = []
    mutations = mutations_str.split(",")
    for mutation in mutations:
        ref_base = mutation[0]
        alt_base = mutation[-1]
        position = int(mutation[1:-1])
        nuc_mutation_objects.append(
            NucMutation(
                strain=strain,
                position=position,
                reference_base=ref_base,
                altered_base=alt_base,
            )
        )
    return nuc_mutation_objects
