<template>
	<section class="overflow-hidden">
		<div class="has-background-design z-0 overflow-clip pt-10">
			<Tabs value="0" class="mx-10 max-w-[70vw] min-w-[70vw]">
				<TabList pt:tabList="!rounded-t-md" pt:activeBar="!bg-white h-0.5">
					<Tab value="0" class="flex items-center gap-2 grow" pt:root="text-white">
						<i class="pi pi-plus-circle pt-1" />
						<span class="font-bold whitespace-nowrap">Metadata Requirement</span>
					</Tab>
					<Tab value="1" class="flex items-center gap-2 grow" pt:root="text-white">
						<i class="pi pi-plus-circle pt-1" />
						<span class="font-bold whitespace-nowrap">Sequence Requirement</span>
					</Tab>
				</TabList>

				<TabPanels>
					<TabPanel value="0" class="m-0 max-h-96 overflow-y-scroll">
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
					</TabPanel>

					<TabPanel value="1" class="m-0 max-h-96 overflow-y-scroll">
						<div class="mb-20">
							<div class="flex justify-center my-4">
								<h3 class="text-lg font-semibold mr-3">Supported file formats:</h3>
								<span
									class="bg-blue-100 text-blue-800 text-sm font-medium mr-2 px-2 pt-1 rounded dark:bg-blue-900 dark:text-blue-300"
								>
									fa
								</span>
								<span
									class="bg-blue-100 text-blue-800 text-sm font-medium mr-2 px-2 pt-1 rounded dark:bg-blue-900 dark:text-blue-300"
								>
									fasta
								</span>
								<span
									class="bg-blue-100 text-blue-800 text-sm font-medium mr-2 px-2 pt-1 rounded dark:bg-blue-900 dark:text-blue-300"
								>
									fna
								</span>
							</div>

							<div class="flex justify-center">
								<h3 class="text-xl font-semibold mb-2">
									Mandatory requirements for the sequence are mentioned below:
								</h3>
							</div>

							<div class="flex flex-col items-center">
								<h6>
									- Metadata-Sequence Header Matching: The strain column in your metadata file
									must exactly match the corresponding sequence headers in the FASTA file. Any
									discrepancies may result in errors or incomplete processing of your data.
								</h6>
								<h6>
									- The sequence data must be submitted as a multi-FASTA file. Ensure that all
									sequence headers in the FASTA file are unique and correspond directly to the
									strain names listed in the metadata. Each sequence should be represented in the
									following format:
								</h6>
							</div>

							<div class="pt-2 pb-7 pl-1">
								>sequence1<br />
								NTAAAGGTTTATACCTTCCCAGGTAACAAACC......<br />
								>sequence2<br />
								NTAAAGGTTTATACCTTCCCAGGTAACAAACC......<br />
								>sequence3<br />
								NTAAAGGTTTATACCTTCCCAGGTAACAAACC......<br />
								...<br />
								>sequence300<br />
								NTAAAGGTTTATACCTTCCCAGGTAACAAACC......<br />
							</div>
						</div>
					</TabPanel>
				</TabPanels>
			</Tabs>
		</div>

		<div
			class="shadow-md rounded-md bg-white dark:bg-[#393939] z-40 -mt-24 relative pt-4 mx-4 md:mx-12 lg:mx-24 xl:mx-48"
		>
			<div class="mx-4 grid lg:grid-cols-2 md:grid-cols-2 md:gap-4 lg:gap-4 grid-cols-1">
				<UploadMetadata @verification_status="verifyMetadata" />
				<UploadSequence @verification_status="verifySequence" />
			</div>

			<div class="flex justify-center gap-12 py-4">
				<Button
					raised
					rounded
					class="!px-10"
					severity="warn"
					@click="download_data"
					v-if="show_qc_check_result"
					:label="`Download QC Failed Data (${uniq_errors.length})`"
				/>
				<Button
					raised
					rounded
					class="!px-10"
					severity="success"
					@click="RunChecks"
					label="Verify Data"
					:disabled="enable_verify"
					v-if="!show_qc_check_result"
				/>
				<Button
					raised
					rounded
					class="!px-10"
					severity="success"
					v-if="show_qc_check_result"
					:label="`Upload QC Passed Data (${cleared_data.length})`"
				/>
			</div>
		</div>

		<div class="mx-8 mt-16 mb-8" v-if="!enable_verify && show_qc_check_result">
			<Accordion :value="accordion_open_index" expandIcon="pi pi-chevron-down" multiple>
				<AccordionPanel
					:key="index"
					:value="index"
					:disabled="!qc_check.data.length"
					v-for="(qc_check, index) in all_qc_checks"
					:style="{ animationDelay: `${index * 250}ms` }"
					class="opacity-0 translate-y-4 transition-all duration-500 ease-in-out animate-fadein"
				>
					<AccordionHeader>
						<span class="flex items-center gap-2 w-full">
							<span class="font-bold whitespace-nowrap">{{ qc_check.name }}</span>
							<div class="ml-auto mr-2">
								<Badge
									severity="warn"
									class="ml-auto mr-2"
									:value="qc_check.error"
									v-if="!qc_check.verification"
								/>
								<span
									v-if="qc_check.verification"
									class="bg-green-100 text-green-800 text-xs font-medium px-2 py-1 rounded-full dark:bg-green-900 dark:text-green-300"
								>
									<i class="pi pi-check text-xs" /> Passed
								</span>

								<span
									v-else
									class="bg-red-100 text-red-800 text-xs font-medium px-2 py-1 rounded-full dark:bg-red-900 dark:text-red-300"
								>
									<i class="pi pi-times text-xs" /> Failed
								</span>
							</div>
						</span>
					</AccordionHeader>

					<AccordionContent v-if="qc_check.data.length">
						<ul
							class="w-full text-sm font-medium text-gray-900 bg-white border border-gray-200 rounded-lg dark:bg-gray-700 dark:border-gray-600 dark:text-white"
						>
							<li
								:key="index"
								v-for="(item, index) in qc_check.data"
								:class="{
									'rounded-t-lg': index == 0,
									'rounded-b-lg border-b-0': index == qc_check.data.length - 1,
								}"
								class="px-4 py-2 border-b border-gray-200 dark:border-gray-600 hover:bg-gray-100 hover:text-gray-950"
							>
								{{ item }}
							</li>
						</ul>
					</AccordionContent>
				</AccordionPanel>
			</Accordion>
		</div>
	</section>
</template>

<script setup>
import JSZip from 'jszip'
const { default: Fasta } = await import('biojs-io-fasta')
import { json2csv } from 'json-2-csv'
import { forEach, groupBy, filter, map, difference, keys, flatten, uniq } from 'lodash'

const dayjs = useDayjs()
const metadata = ref(null)
const sequence = ref(null)
const show_qc_check_result = ref(false)
const accordion_open_index = computed(() => map(all_qc_checks.value, (d, i) => i))
const enable_verify = computed(() => {
	const metadataLength = metadata.value ? metadata.value.length : 0
	const sequenceLength = sequence.value ? keys(sequence.value).length : 0

	const return_value = !(metadataLength > 0 && sequenceLength > 0)

	if (return_value) {
		show_qc_check_result.value = false
		all_qc_checks.value = all_qc_checks.value.slice(0, 4)
	}

	return return_value
})
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
const uniq_errors = ref([])
const cleared_data = ref([])
const all_qc_checks = ref([
	{
		data: [],
		error: null,
		show: false,
		verification: false,
		name: 'Metadata file format check',
	},
	{
		data: [],
		error: null,
		show: false,
		verification: false,
		name: 'Metadata file structure check',
	},
	{
		data: [],
		error: null,
		show: false,
		verification: false,
		name: 'Sequence file format check',
	},
	{
		data: [],
		error: null,
		show: false,
		verification: false,
		name: 'Sequence file structure check',
	},
])

const verifyMetadata = (data) => {
	all_qc_checks.value[0].verification = data.verification
	all_qc_checks.value[1].verification = data.verification
	metadata.value = data.data
}

const verifySequence = (data) => {
	all_qc_checks.value[2].verification = data.verification
	all_qc_checks.value[3].verification = data.verification
	sequence.value = data.data
}

const RunChecks = () => {
	// Metadata based checks
	check_collection_date()
	find_duplicate_metadata()
	// Sequence based checks
	find_duplicate_sequence()
	match_metadata_with_sequence()

	show_qc_check_result.value = true
	uniq_errors.value = uniq(flatten(map(all_qc_checks.value, (d) => d.data)))
	cleared_data.value = difference(uniq(map(metadata.value, (d) => d['Virus name'])), uniq_errors.value)
	console.log('🚀 ~ RunChecks ~ uniq_errors:', uniq_errors, cleared_data)
}

const check_collection_date = () => {
	const collection_date_error = []
	const collection_date_early = []
	const empty_collection_date = []
	const collection_date_future = []

	forEach(metadata.value, (d) => {
		// First Metadata Verification - Empty Collection dates
		if (d['Collection date'] === '') {
			empty_collection_date.push(d['Virus name'])
		} else {
			const collectionDate = dayjs(
				d['Collection date'],
				['DD-MM-YYYY', 'DD/MM/YYYY', 'YYYY-MM-DD', 'YYYY/MM/DD'],
				true,
			)

			// Second Metadata Verification - Error in collection dates
			if (!collectionDate.isValid()) {
				collection_date_error.push(d['Virus name'])
			}
			// Third Metadata Verification - Date cannot be before December 31, 2019 as that was the first collection of the sample in Wuhan
			if (collectionDate.isBefore('2019-12-31', 'date')) {
				collection_date_early.push(d['Virus name'])
			}
			// Fourth Metadata Verification - Date cannot be in the future
			if (collectionDate.isAfter(dayjs(), 'day')) {
				collection_date_future.push(d['Virus name'])
			}
		}
	})
	all_qc_checks.value.push({
		show: false,
		data: empty_collection_date,
		error: empty_collection_date.length,
		name: 'Empty collection date for following',
		verification: !empty_collection_date.length,
	})
	all_qc_checks.value.push({
		show: false,
		data: collection_date_error,
		error: collection_date_error.length,
		verification: !collection_date_error.length,
		name: 'Error in collection date format (DD-MM-YYYY, DD/MM/YYYY, YYYY-MM-DD, YYYY/MM/DD)',
	})
	all_qc_checks.value.push({
		show: false,
		data: collection_date_early,
		error: collection_date_early.length,
		verification: !collection_date_early.length,
		name: 'Collection date for following earlier than 2019-12-31',
	})
	all_qc_checks.value.push({
		show: false,
		data: collection_date_future,
		error: collection_date_future.length,
		name: 'Collection date in the future',
		verification: !collection_date_future.length,
	})
}

const find_duplicate_metadata = () => {
	// Fifth Metadata Verification - Duplicated metadata which may or may not contain different values of other fields
	const virusNameCounts = groupBy(metadata.value, (d) => d['Virus name'])
	const duplicateVirusNames = map(
		filter(virusNameCounts, (group) => group.length > 1),
		(group) => group[0]['Virus name'],
	)

	all_qc_checks.value.push({
		show: false,
		data: duplicateVirusNames,
		name: 'Duplicate metadata',
		error: duplicateVirusNames.length,
		verification: !duplicateVirusNames.length,
	})
}

const find_duplicate_sequence = () => {
	// First Sequence Verification - Duplicated sequence which may or may not contain different fasta
	const virusNameCounts = groupBy(sequence.value, (d) => d['name'])
	const duplicateVirusNames = map(
		filter(virusNameCounts, (group) => group.length > 1),
		(group) => group[0]['Virus name'],
	)

	all_qc_checks.value.push({
		show: false,
		data: duplicateVirusNames,
		name: 'Duplicate sequence',
		error: duplicateVirusNames.length,
		verification: !duplicateVirusNames.length,
	})
}

const match_metadata_with_sequence = () => {
	const metadata_virus_name = map(metadata.value, (d) => d['Virus name'])
	const sequence_virus_name = map(sequence.value, (d) => d['name'])

	// Sixth Metadata Verification - Matching all metadata with the sequence, so for every metadata there must be a sequence
	const missing_sequences = difference(metadata_virus_name, sequence_virus_name)
	all_qc_checks.value.push({
		show: false,
		data: missing_sequences,
		name: 'Sequence Missing',
		error: missing_sequences.length,
		verification: !missing_sequences.length,
	})
	// Second Sequence Verification - Matching all sequence with the metadata, so for every sequence there must be a metadata
	const missing_metadata = difference(sequence_virus_name, metadata_virus_name)
	all_qc_checks.value.push({
		show: false,
		data: missing_metadata,
		name: 'Metadata Missing',
		error: missing_metadata.length,
		verification: !missing_metadata.length,
	})
}

const download_blob = (blob, filename) => {
	const link = document.createElement('a')
	link.href = URL.createObjectURL(blob)
	link.download = filename
	document.body.appendChild(link)
	link.click()
	document.body.removeChild(link)
}

const download_data = async () => {
	try {
		const zip = new JSZip()

		const qcFailedData = filter(metadata.value, (d) => uniq_errors.value.includes(d['Virus name']))
		const qcFailedSequences = filter(sequence.value, (d) => uniq_errors.value.includes(d['name']))

		// Add files to ZIP
		zip.file('qc_failed_metadata.tsv', json2csv(qcFailedData, { delimiter: { field: '\t' } }))
		zip.file('qc_failed_sequences.fasta', Fasta.write(qcFailedSequences))

		// Show loader
		// loader.value = true;

		// Generate ZIP and download
		const blob = await zip.generateAsync({
			type: 'blob',
			compression: 'DEFLATE',
			compressionOptions: { level: 5 },
		})

		download_blob(blob, 'qc_failed_metadata.zip')
	} catch (error) {
		console.error('Error downloading data:', error)
	}
}
</script>

<style>
.has-background-design {
	width: 120vw;
	display: flex;
	min-height: 60vh;
	margin-left: -10vw;
	position: relative;
	align-items: center;
	flex-direction: column;
	justify-content: flex-end;
	border-radius: 0 0 50% 50%;
	/*background-image: linear-gradient(0deg,#bfeaff,#80caff);*/
	background-image: linear-gradient(0deg, rgb(0, 60, 89), rgb(0, 74, 127));
}
</style>
