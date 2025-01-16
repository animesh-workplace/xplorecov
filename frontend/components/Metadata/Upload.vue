<template>
	<section>
		<Skeleton v-if="isLoading" width="100%" height="5rem" class="mb-2" />
		<div :class="{ hidden: isLoading }">
			<client-only>
				<file-pond
					max-files="1"
					credits="false"
					@init="RemoveLoader"
					@addfile="HandleFile"
					class="cursor-pointer"
					ref="metadata_filepond"
					:allow-multiple="false"
					@removefile="RemoveFile"
					v-model:files="metadata_file"
					label-idle="
                        <span class='cursor-pointer font-bold'>
                            Drag & Drop your Metadata or
                        </span>
                        <span class='cursor-pointer underline font-bold'>
                            Browse
                        </span>
                    "
				/>
			</client-only>
		</div>
	</section>
</template>

<script setup>
/*
This component provides a file upload interface specifically designed for handling metadata files related to SARS-CoV-2 sequences.
It uses the Vue FilePond library for an enhanced file upload experience, supporting drag-and-drop functionality and providing visual feedback.

### Key Features:
1. **File Upload with Validation:**
   - Allows only a single file upload (`max-files="1"`, `allow-multiple="false`).
   - Restricts file types to CSV, TSV, or text files.
   - Validates the presence of required metadata columns in the uploaded file.

2. **Real-time Feedback:**
   - Shows detailed error messages if the file type is invalid or required columns are missing.

3. **File Processing Logic:**
   - Converts the uploaded file to a cleaned JSON format (removes carriage returns and empty rows).
   - Verifies the metadata structure by checking for the presence of all required columns.
   - Updates a `metadata_verification` state to track the results of various validation checks.

4. **Error Handling:**
   - Leverages the `usePondFileError` composable to provide user-friendly error notifications for issues like invalid file types or missing columns.

### Dependencies:
- **Vue FilePond Plugins:**
  - `FilePondPluginFileValidateType` for file type validation.
  - `FilePondPluginImagePreview` for potential image preview capabilities (not actively used here but included).
- **Lodash Functions:** Used for cleaning JSON keys, filtering data, and checking for missing columns.
- **CSV Conversion Library (`csvjson-csv2json`):** Converts CSV/TSV files into JSON for further processing.

### Event Handlers:
- `RemoveLoader`: Disables the loading state once the component initializes.
- `HandleFile`: Processes the uploaded file, performs validations, and notifies the user of errors or success.
- `RemoveFile`: Resets all verification checks when the file is removed.
*/

import csv2json from 'csvjson-csv2json'
import vueFilePond from 'vue-filepond'
import { usePondFileError } from '@/composables/usePondFileError'
import { map, cloneDeepWith, isString, difference, join, filter } from 'lodash-es'
import FilePondPluginImagePreview from 'filepond-plugin-image-preview/dist/filepond-plugin-image-preview.esm.js'
import FilePondPluginFileValidateType from 'filepond-plugin-file-validate-type/dist/filepond-plugin-file-validate-type.esm.js'

import 'filepond/dist/filepond.min.css'
import 'filepond-plugin-image-preview/dist/filepond-plugin-image-preview.min.css'

const FilePond = vueFilePond(FilePondPluginFileValidateType, FilePondPluginImagePreview)

// Reactive state for files
const isLoading = ref(true)
const metadata_file = ref(null)
const metadata_filepond = ref(null)
const emit = defineEmits(['verification_status'])
const { handlePondFileError } = usePondFileError()

// Remove Loader when the page is loaded
const RemoveLoader = () => {
	isLoading.value = false
}

// Removes \r from text
const cleanJSON = (array) => {
	return map(array, (obj) =>
		cloneDeepWith(obj, (value) => {
			if (isString(value)) {
				return value.replace(/\r/g, '')
			}
		}),
	)
}

// Remove the rows which are empty that donot contain any Virus name
const removeEmptyRows = (obj) => {
	return filter(obj, (d) => d['Virus name'] != '')
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
	fileReader.onload = (e) => {
		const metadata = e.target.result
		// Converting the file to cleaned json
		metadata_json = removeEmptyRows(cleanJSON(csv2json(metadata, { parseNumbers: true })))

		// Second Check:  Validation of the columns whether all necessary columns are present or not
		const checking_columns_result = checkColumns(metadata_json)
		if (checking_columns_result.return_value) {
			emit('verification_status', {
				verification: true,
				data: metadata_json,
			})
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
	emit('verification_status', {
		verification: false,
		data: null,
	})
}
</script>

<style></style>
