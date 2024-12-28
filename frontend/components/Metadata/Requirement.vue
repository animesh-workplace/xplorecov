<template>
	<div>
		<div class="flex justify-center my-4">
			<h3 class="text-lg font-semibold mr-3">Supported file formats:</h3>

			<span
				class="bg-blue-100 text-blue-800 text-sm font-medium mr-2 px-2 pt-1 rounded dark:bg-blue-900 dark:text-blue-300"
			>
				tsv
			</span>
			<span
				class="bg-blue-100 text-blue-800 text-sm font-medium mr-2 px-2 pt-1 rounded dark:bg-blue-900 dark:text-blue-300"
			>
				csv
			</span>
			<span
				class="bg-blue-100 text-blue-800 text-sm font-medium mr-2 px-2 pt-1 rounded dark:bg-blue-900 dark:text-blue-300"
			>
				txt
			</span>
		</div>

		<div class="flex justify-center">
			<h3 class="text-xl font-semibold mb-2">
				Mandatory/Optional headers for the metadata are mentioned below:
			</h3>
		</div>

		<DataTable :value="metadata_requirements" class="mb-32 mr-4">
			<Column field="name" header="Header" />
			<Column field="required" header="Mandatory">
				<template #body="slotProps">
					<span
						v-if="slotProps.data.required"
						class="bg-red-100 text-red-800 text-xs font-medium px-2 py-1 rounded-full dark:bg-red-900 dark:text-red-300"
					>
						mandatory
					</span>
					<span
						v-else
						class="bg-gray-100 text-gray-800 text-xs font-medium px-2 py-1 rounded-full dark:bg-gray-700 dark:text-gray-300"
					>
						optional
					</span>
				</template>
			</Column>
			<Column field="description" header="Description" />
		</DataTable>
	</div>
</template>

<script setup>
const metadata_requirements = ref([
	{
		name: 'Virus name',
		required: true,
		description: 'e.g. hCoV-19/India/1234/2021 (Must be FASTA-Header from the FASTA file)',
	},
	{
		name: 'Type',
		required: true,
		description: 'default must remain betacoronavirus',
	},
	{
		name: 'Passage details/history',
		required: true,
		description: 'e.g. Original, Vero',
	},
	{
		name: 'Collection date',
		required: true,
		description: 'Format [ DD-MM-YYYY, DD/MM/YYYY, YYYY-MM-DD, YYYY/MM/DD ]',
	},
	{
		name: 'Continent',
		required: true,
		description: 'Sample collected from which continent',
	},
	{
		name: 'Country',
		required: true,
		description: 'Sample collected from which country',
	},
	{
		name: 'State',
		required: true,
		description: 'Sample collected from which state',
	},
	{
		name: 'District',
		required: true,
		description: 'Sample collected from which district',
	},
	{
		name: 'Additional location information',
		required: false,
		description: 'e.g. Cruise Ship, Convention, Live animal market',
	},
	{
		name: 'Host',
		required: true,
		description: 'e.g. Human, Environment, Canine, Manis javanica, Rhinolophus affinis, etc',
	},
	{
		name: 'Additional host information',
		required: false,
		description: 'e.g. Patient infected while traveling in …',
	},
	{
		name: 'Sampling Strategy',
		required: false,
		description: 'e.g. Sentinel surveillance (ILI, ARI, SARI), Non-sentinel-surveillance',
	},
	{
		name: 'Gender',
		required: true,
		description: 'Male, Female, or unknown',
	},
	{
		name: 'Patient age',
		required: true,
		description: 'e.g. 65 or 7 months, or unknown',
	},
	{
		name: 'Patient status',
		required: true,
		description: 'e.g. Hospitalized, Released, Live, Deceased, or unknown',
	},
	{
		name: 'Specimen source',
		required: false,
		description: 'e.g. Sputum, Oro-pharyngeal swab, Blood, Tracheal swab, Other',
	},
	{
		name: 'Outbreak',
		required: false,
		description: 'Date, Location e.g. type of gathering, Family cluster, etc.',
	},
	{
		name: 'Last vaccinated',
		required: false,
		description: 'provide details if applicable',
	},
	{
		name: 'Treatment',
		required: false,
		description: 'Include drug name, dosage',
	},
	{
		name: 'Sequencing technology',
		required: true,
		description: 'e.g.	Illumina Miseq, Sanger, Nanopore MinION, Ion Torrent, etc.',
	},
	{
		name: 'Assembly method',
		required: false,
		description: 'e.g. CLC Genomics Workbench 12, SPAdes/MEGAHIT v1.2.9, etc.',
	},
	{
		name: 'Coverage',
		required: false,
		description: 'e.g. 70x, 1000x, 10000x',
	},
	{
		name: 'Originating lab',
		required: true,
		description: 'Where the clinical specimen or virus isolate was first obtained',
	},
	{
		name: 'Originating lab address',
		required: true,
		description: '',
	},
	{
		name: 'Sample ID given by the originating lab',
		required: false,
		description: '',
	},
	{
		name: 'Submitting lab',
		required: true,
		description: 'Where sequence data have been generated and submitted to GISAID',
	},
	{
		name: 'Submitting lab address',
		required: true,
		description: '',
	},
	{
		name: 'Sample ID given by the submitting lab',
		required: false,
		description: '',
	},
	{
		name: 'Authors',
		required: true,
		description: 'a comma separated list of Authors with First followed by Last Name',
	},
])
</script>

<style scoped></style>
