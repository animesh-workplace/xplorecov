<template>
	<section>
		<Skeleton v-if="isLoading" width="100%" height="5rem" class="mb-2" />
		<div :class="{ hidden: isLoading }">
			<client-only>
				<file-pond
					max-files="1"
					credits="false"
					:allow-multiple="false"
					ref="sequence_filepond"
					@init="RemoveLoader"
					@addfile="HandleFile"
					@removefile="RemoveFile"
					v-model:files="sequence_file"
					label-idle="
                        <span class='is-family-primary has-text-weight-semibold has-text-grey-dark is-clickable'>
                            Drag & Drop your Multi-Sequence fasta or
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
import vueFilePond from 'vue-filepond'
import { forEach, mapValues, mapKeys } from 'lodash'
import { usePondFileError } from '@/composables/usePondFileError'
import FilePondPluginFileValidateType from 'filepond-plugin-file-validate-type/dist/filepond-plugin-file-validate-type.esm.js'

import 'filepond/dist/filepond.min.css'

const FilePond = vueFilePond(FilePondPluginFileValidateType)

// Reactive state for files
const isLoading = ref(true)
const sequence_file = ref(null)
const sequence_filepond = ref(null)
const { handlePondFileError } = usePondFileError()
const sequence_verification = ref([
	{ name: 'Sequence Fasta file format check', verification: false },
	{ name: 'Sequence Fasta file structure check', verification: false },
	{ name: 'Sequence Fasta file metadata check', verification: false },
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
			sequence_verification.value[1].verification = true
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
		sequence_verification.value[0].verification = true
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
	forEach(sequence_verification.value, (item) => (item.verification = false))
}
</script>

<style></style>
