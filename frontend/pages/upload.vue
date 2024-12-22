<template>
	<section>
		<div class="mx-8 mt-8 mb-4 grid lg:grid-cols-2 lg:gap-14 grid-cols-1">
			<UploadMetadata @verification_status="verifyMetadata" />
			<UploadSequence @verification_status="verifySequence" />
		</div>

		<div class="flex justify-center">
			<Button
				raised
				rounded
				class="!px-10"
				severity="success"
				@click="RunChecks"
				label="Verify Data"
				:disabled="enable_verify"
			/>
		</div>

		<div class="mx-8 my-8" v-if="!enable_verify && show_qc_check_result">
			<Accordion :value="accordion_open_index" expandIcon="pi pi-chevron-down" multiple>
				<AccordionPanel
					:value="index"
					v-for="(qc_check, index) in all_qc_checks"
					:key="index"
					:disabled="!qc_check.data.length"
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
import { forEach, groupBy, filter, map, difference, keys } from 'lodash'
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
	console.log('🚀 ~ missing_sequences:', missing_sequences)
	all_qc_checks.value.push({
		show: false,
		data: missing_sequences,
		name: 'Sequence Missing',
		error: missing_sequences.length,
		verification: !missing_sequences.length,
	})
	// Second Sequence Verification - Matching all sequence with the metadata, so for every sequence there must be a metadata
	const missing_metadata = difference(sequence_virus_name, metadata_virus_name)
	console.log('🚀 ~ missing_metadata:', missing_metadata)
	all_qc_checks.value.push({
		show: false,
		data: missing_metadata,
		name: 'Metadata Missing',
		error: missing_metadata.length,
		verification: !missing_metadata.length,
	})
}
</script>
