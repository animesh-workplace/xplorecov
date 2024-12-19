<template>
	<section>
		<Skeleton v-if="isLoading" width="100%" height="5rem" class="mb-2" />
		<div :class="{ hidden: isLoading }">
			<client-only>
				<file-pond
					max-files="1"
					credits="false"
					:allow-multiple="false"
					ref="metadata_filepond"
					@init="RemoveLoader"
					@addfile="HandleFile"
					@removefile="RemoveFile"
					v-model:files="metadata_file"
					label-idle="
                        <span class='is-family-primary has-text-weight-semibold has-text-grey-dark is-clickable'>
                            Drag & Drop your Metadata or
                        </span>
                        <span class='is-family-primary has-text-weight-semibold has-text-grey-dark is-clickable'>
                            Browse
                        </span>
                    "
				/>
			</client-only>
		</div>
	</section>
</template>

<script setup>
import csv2json from 'csvjson-csv2json'
import vueFilePond from 'vue-filepond'
import { usePondFileError } from '@/composables/usePondFileError'
import { forEach, mapValues, mapKeys, difference, join } from 'lodash'
import FilePondPluginImagePreview from 'filepond-plugin-image-preview/dist/filepond-plugin-image-preview.esm.js'
import FilePondPluginFileValidateType from 'filepond-plugin-file-validate-type/dist/filepond-plugin-file-validate-type.esm.js'

import 'filepond/dist/filepond.min.css'
import 'filepond-plugin-image-preview/dist/filepond-plugin-image-preview.min.css'

const FilePond = vueFilePond(FilePondPluginFileValidateType, FilePondPluginImagePreview)

// Reactive state for files
const isLoading = ref(true)
const metadata_file = ref(null)
const metadata_filepond = ref(null)
const { handlePondFileError } = usePondFileError()
const metadata_verification = ref([
	{ name: 'Metadata format check', verification: false },
	{ name: 'Metadata structure check', verification: false },
	{ name: 'Metadata sequence check', verification: false },
])

// Methods
const RemoveLoader = () => {
	isLoading.value = false
}

// Removes \r from text
const cleanJSON = (obj) => {
	const cleanedKeys = mapKeys(obj, (value, key) => key.replace(/\r/g, '')) // Clean keys
	return mapValues(cleanedKeys, (value) => (typeof value === 'string' ? value.replace(/\r/g, '') : value))
}

// Checks if the columns contain required columns or not
const checkColumns = (metadataJson) => {
	const keys = Object.keys(metadataJson[0]) // Get keys from the first object
	const requiredColumns = [
		'Virus name',
		'Type',
		'Passage details/history',
		'Collection date',
		'Country',
		'State',
		'District',
		'Location',
		'Additional location information',
		'Host',
		'Additional host information',
		'Gender',
		'Patient age',
		'Patient status',
		'Specimen source',
		'Outbreak',
		'Last vaccinated',
		'Treatment',
		'Sequencing technology',
		'Assembly method',
		'Coverage',
		'Originating lab',
		'Originating lab address',
		'Submitting lab',
		'Submitting lab address',
		'Sample ID given by the submitting lab',
		'Authors',
	]

	// Calculate missing columns
	const missingColumns = difference(requiredColumns, keys)

	// Return results
	return {
		missing: missingColumns,
		return_value: missingColumns.length === 0 ? 1 : 0,
	}
}

const HandleFile = (error, file) => {
	let metadata_json
	const file_type_accepted = ['text', 'csv', 'tsv', 'txt']
	const fileReader = new FileReader()
	fileReader.onload = function (e) {
		const metadata = e.target.result
		// Converting the file to cleaned json
		metadata_json = cleanJSON(csv2json(metadata, { parseNumbers: true }))

		// Second Check:  Validation of the columns whether all necessary columns are present or not
		const checking_columns_result = checkColumns(metadata_json)
		if (checking_columns_result.return_value) {
			metadata_verification.value[1].verification = true
		} else {
			const missing_text = join(checking_columns_result.missing, '\n')
			push.error({
				title: 'Following columns are missing:',
				message: missing_text,
				duration: undefined,
			})
			handlePondFileError(
				metadata_filepond.value,
				'processing-warn',
				'Columns not matching',
				'Check metadata requirements',
			)
		}
	}

	// First Check: Validation of the file extension and providing proper notification
	if (file_type_accepted.includes(file.fileExtension)) {
		metadata_verification.value[0].verification = true
		fileReader.readAsText(file.file)
	} else {
		push.error('File is of invalid type, expects csv, tsv or text/txt')
		handlePondFileError(
			metadata_filepond.value,
			'processing-error',
			'File is of invalid type',
			'Expects csv, tsv, txt file',
		)
	}
}

const RemoveFile = () => {
	forEach(metadata_verification.value, (item) => (item.verification = false))
}
</script>

<style>
[data-filepond-item-state='processing-warn'] .filepond--item-panel {
	background-color: #d97a37;
}

[data-filepond-item-state*='processing-warn'] .filepond--panel,
[data-filepond-item-state*='processing-warn'] .filepond--file-wrapper {
	-webkit-animation: shake 0.65s linear both;
	animation: shake 0.65s linear both;
}
</style>
