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
					:allow-multiple="false"
					ref="sequence_filepond"
					@removefile="RemoveFile"
					v-model:files="sequence_file"
					label-idle="
						<span class='cursor-pointer font-bold'>
                            Drag & Drop your Multi-Sequence fasta or
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
This component provides a file upload interface tailored for multi-sequence FASTA files.
It leverages the Vue FilePond library to enable a drag-and-drop file upload experience with real-time feedback and validation checks.

### Key Features:
1. **File Upload with Validation:**
   - Allows only a single file upload (`max-files="1"`, `allow-multiple="false`).
   - Restricts file types to FASTA files (`fasta`, `fa`, `fna`).
   - Performs multiple verification checks, including:
     - File format validation.
     - File structure verification.

2. **Real-time Feedback:**
   - Provides detailed error messages if the file type is invalid or if structural issues are detected in the FASTA file.

3. **File Processing Logic:**
   - Parses and cleans FASTA files using the `biojs-io-fasta` library.
   - Converts uploaded sequences into a cleaned JSON format by removing unwanted carriage returns (`\r`).
   - Updates the `sequence_verification` state to track the results of file verification steps.

4. **Error Handling:**
   - Integrates the `usePondFileError` composable for user-friendly error notifications.
   - Catches issues related to unsupported file types or invalid FASTA file structures and notifies the user with clear messages.

### Dependencies:
- **Vue FilePond Plugins:**
  - `FilePondPluginFileValidateType` for file type validation.
- **Lodash Functions:** Utilized for cleaning JSON keys and resetting validation states.
- **biojs-io-fasta Library:** Parses multi-sequence FASTA files into a usable JSON format for validation and processing.

### Event Handlers:
- `RemoveLoader`: Disables the loading state once the component initializes.
- `HandleFile`: Processes the uploaded file, validates the file type and structure, and provides error notifications for invalid files.
- `RemoveFile`: Resets all verification checks when a file is removed.
*/

import vueFilePond from 'vue-filepond'
import { map, cloneDeepWith, isString } from 'lodash'
import { usePondFileError } from '@/composables/usePondFileError'
import FilePondPluginFileValidateType from 'filepond-plugin-file-validate-type/dist/filepond-plugin-file-validate-type.esm.js'

import 'filepond/dist/filepond.min.css'

const FilePond = vueFilePond(FilePondPluginFileValidateType)

// Reactive state for files
const isLoading = ref(true)
const sequence_file = ref(null)
const sequence_filepond = ref(null)
const emit = defineEmits(['verification_status'])
const { handlePondFileError } = usePondFileError()

// Methods
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

const HandleFile = (error, file) => {
	let sequence_json
	const file_type_accepted = ['fasta', 'fa', 'fna']
	const fileReader = new FileReader()
	fileReader.onload = async (e) => {
		const sequence = e.target.result
		// Second Check: Whether converting the file to cleaned json is possible or not due to file structure issues
		const { default: Fasta } = await import('biojs-io-fasta')
		try {
			sequence_json = cleanJSON(Fasta.parse(sequence))
			emit('verification_status', {
				verification: true,
				data: sequence_json,
			})
		} catch (err) {
			console.log(err)
			push.error('Invalid fasta file')
			handlePondFileError(
				sequence_filepond.value,
				'processing-warn',
				'Invalid fasta file',
				'Unable to parse the fasta file',
			)
		}
	}

	// First Check: Validation of the file extension and providing proper notification
	if (file_type_accepted.includes(file.fileExtension)) {
		fileReader.readAsText(file.file)
	} else {
		push.error('File is of invalid type, expects fasta or fa or fna')
		handlePondFileError(
			sequence_filepond.value,
			'processing-error',
			'File is of invalid type',
			'Expects fasta or fa or fna',
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
