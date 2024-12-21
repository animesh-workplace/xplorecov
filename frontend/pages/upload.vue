<template>
	<section>
		<div class="mx-8 mt-8 mb-4 grid lg:grid-cols-2 lg:gap-14 grid-cols-1">
			<UploadMetadata @verification_status="verifyMetadata" />
			<UploadSequence @verification_status="verifySequence" />
		</div>

		<div class="flex justify-center">
			<Button raised rounded class="!px-10" severity="success" label="Verify Data" @click="RunChecks" />
		</div>
	</section>
</template>

<script setup>
import { forEach, groupBy, filter, map, difference } from 'lodash'
const dayjs = useDayjs()
const metadata = ref([])
const sequence = ref([])

const all_qc_checks = ref([
	{
		name: 'Metadata file format check',
		verification: false,
		data: [],
		error: null,
		show: false,
	},
	{
		name: 'Metadata file structure check',
		verification: false,
		data: [],
		error: null,
		show: false,
	},
	{
		name: 'Sequence file format check',
		verification: false,
		data: [],
		error: null,
		show: false,
	},
	{
		name: 'Sequence file structure check',
		verification: false,
		data: [],
		error: null,
		show: false,
	},
	// { name: 'Missing Metadata/Sequence check', verification: false, data: [], error: null, show: false },
	// { name: 'Duplicate check', verification: false, data: [], error: null, show: false },
	// { name: 'Already present check', verification: false, data: [], error: null, show: false },
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
</script>
