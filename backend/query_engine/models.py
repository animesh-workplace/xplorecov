from django.db import models
from django.contrib.postgres.fields import ArrayField
from django.core.validators import MinValueValidator, MaxValueValidator

# Create your models here.
# models.py

class VirusStrain(models.Model):
    PANGOLIN_QC_CHOICES = [
        ('fail', 'Fail'),
        ('pass', 'Pass'),
    ]
    NEXTCLADE_QC_CHOICES = [
        ('good', 'Good'),
        ('mediocre', 'Mediocre'),
        ('bad', 'Bad'),
    ]
    strain = models.CharField(max_length=255)
    gisaid_epi_isl = models.CharField(max_length=50, unique=True)
    genbank_accession = models.CharField(max_length=50, null=True, blank=True)
    collection_date = models.DateField(null=True, blank=True)
    host = models.CharField(max_length=100, null=True, blank=True)
    age = models.IntegerField(null=True, blank=True)
    sex = models.CharField(max_length=10, null=True, blank=True)
    coverage = models.DecimalField(
        max_digits=10, decimal_places=5, null=True, blank=True,
        validators=[MinValueValidator(0), MaxValueValidator(1)]
    )
    nextclade_clade = models.CharField(max_length=50, null=True, blank=True)
    nextclade_lineage = models.CharField(max_length=50, null=True, blank=True)
    pangolin_lineage = models.CharField(max_length=50, null=True, blank=True)
    pangolin_qc_status = models.CharField(max_length=20, null=True, blank=True, choices=PANGOLIN_QC_CHOICES)
    nextclade_qc_status = models.CharField(
        max_length=20, choices=NEXTCLADE_QC_CHOICES, null=True, blank=True
    )    

    class Meta:
        indexes = [
            models.Index(fields=['collection_date']),
            models.Index(fields=['nextclade_clade']),
            models.Index(fields=['pangolin_lineage']),
            models.Index(fields=['nextclade_lineage']),
        ]

    def __str__(self):
        return f"{self.strain} ({self.gisaid_epi_isl})"

class Location(models.Model):
    strain = models.OneToOneField(
        VirusStrain,
        on_delete=models.CASCADE,
        related_name='location'
    )
    continent = models.CharField(max_length=100, null=True, blank=True)
    country = models.CharField(max_length=100, null=True, blank=True)
    state = models.CharField(max_length=100, null=True, blank=True)
    district = models.CharField(max_length=255, null=True, blank=True)
    continent_exposure = models.CharField(max_length=100, null=True, blank=True)
    country_exposure = models.CharField(max_length=100, null=True, blank=True)
    state_exposure = models.CharField(max_length=100, null=True, blank=True)

    class Meta:
        indexes = [
            models.Index(fields=['state']),
            models.Index(fields=['country']),
            models.Index(fields=['district']),
        ]

class AAMutation(models.Model):
    PROTEIN_CHOICES = [
        ('N', 'Nucleocapsid (N)'),
        ('S', 'Spike (S)'),
        ('E', 'Envelope (E)'),
        ('M', 'Membrane (M)'),
        ('ORF1a', 'ORF1a Polyprotein'),
        ('ORF1b', 'ORF1b Polyprotein'),
        ('ORF3a', 'ORF3a Protein'),
        ('ORF3b', 'ORF3b Protein'),
        ('ORF6', 'ORF6 Protein'),
        ('ORF7a', 'ORF7a Protein'),
        ('ORF7b', 'ORF7b Protein'),
        ('ORF8', 'ORF8 Protein'),
        ('ORF9b', 'ORF9b Protein'),
        ('ORF10', 'ORF10 Protein'),
        ('NSP1', 'Non-structural Protein 1 (NSP1)'),
        ('NSP2', 'Non-structural Protein 2 (NSP2)'),
        ('NSP3', 'Non-structural Protein 3 (NSP3)'),
        ('NSP4', 'Non-structural Protein 4 (NSP4)'),
        ('NSP5', 'Main Protease (NSP5)'),
        ('NSP6', 'Non-structural Protein 6 (NSP6)'),
        ('NSP7', 'Non-structural Protein 7 (NSP7)'),
        ('NSP8', 'Non-structural Protein 8 (NSP8)'),
        ('NSP9', 'Non-structural Protein 9 (NSP9)'),
        ('NSP10', 'Non-structural Protein 10 (NSP10)'),
        ('NSP12', 'RNA-dependent RNA Polymerase (NSP12)'),
        ('NSP13', 'Helicase (NSP13)'),
        ('NSP14', '3\'-to-5\' Exonuclease (NSP14)'),
        ('NSP15', 'Endoribonuclease (NSP15)'),
        ('NSP16', '2\'-O-Methyltransferase (NSP16)'),
    ]

    strain = models.ForeignKey(
        VirusStrain,
        on_delete=models.CASCADE,
        related_name='aa_mutations'
    )
    protein = models.CharField(max_length=50, choices=PROTEIN_CHOICES)
    reference_aa = models.CharField(max_length=10)
    position = models.IntegerField()
    altered_aa = models.CharField(max_length=10)

    class Meta:
        indexes = [
            models.Index(fields=['protein', 'position', 'reference_aa', 'altered_aa'])
        ]

class NucMutation(models.Model):
    strain = models.ForeignKey(
        VirusStrain,
        on_delete=models.CASCADE,
        related_name='nuc_mutations'
    )
    reference_base = models.CharField(max_length=1)
    position = models.IntegerField()
    altered_base = models.CharField(max_length=1)

    class Meta:
        indexes = [
            models.Index(fields=['reference_base', 'position', 'altered_base'])
        ]

class StructuralVariation(models.Model):
    strain = models.OneToOneField(
        VirusStrain,
        on_delete=models.CASCADE,
        related_name='structural_variation'
    )
    total_insertions = models.IntegerField(default=0)
    insertion_positions = models.TextField(null=True, blank=True)  # Store as JSON string
    total_frameshifts = models.IntegerField(default=0)
    frameshift_positions = models.TextField(null=True, blank=True)  # Store as JSON string
    total_missing = models.IntegerField(default=0)
    missing_ranges = models.TextField(null=True, blank=True)  # Store as JSON string

    class Meta:
        indexes = [
            models.Index(fields=['total_insertions', 'total_frameshifts', 'total_missing'])
        ]

class QualityMetrics(models.Model):
    strain = models.OneToOneField(
        VirusStrain,
        on_delete=models.CASCADE,
        related_name='quality_metrics'
    )
    total_non_acgtns = models.IntegerField(default=0)
    non_acgtn_positions = models.TextField(null=True, blank=True)  # Store as JSON string
    partially_aliased = models.BooleanField(default=False)
    scorpio_call = models.CharField(max_length=50, null=True, blank=True)
    scorpio_support = models.DecimalField(
        max_digits=5,
        decimal_places=2,
        null=True,
        blank=True
    )
    scorpio_conflict = models.DecimalField(
        max_digits=5,
        decimal_places=2,
        null=True,
        blank=True
    )


# Helper function for mutation search
def parse_aamutation_string(mutation_string):
    """Parse a mutation string into components (e.g., 'N:P13L' -> 
       {'protein': 'N', 'reference': 'P', 'position': 13, 'alternate': 'L'})"""
    try:
        protein, mut = mutation_string.split(':')
        reference = mut[0]  # First character before the position is the reference
        position = int(''.join(filter(str.isdigit, mut)))  # Extract position as an integer
        alternate = mut[-1]  # Last character is the alternate base
        return {
            'protein': protein,
            'reference': reference,
            'position': position,
            'alternate': alternate
        }
    except (ValueError, IndexError):
        return None
