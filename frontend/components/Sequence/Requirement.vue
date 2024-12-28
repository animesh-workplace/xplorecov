<template>
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
				- Metadata-Sequence Header Matching: The strain column in your metadata file must exactly match the
				corresponding sequence headers in the FASTA file. Any discrepancies may result in errors or
				incomplete processing of your data.
			</h6>
			<h6>
				- The sequence data must be submitted as a multi-FASTA file. Ensure that all sequence headers in
				the FASTA file are unique and correspond directly to the strain names listed in the metadata. Each
				sequence should be represented in the following format:
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
</template>

<script setup></script>

<style scoped></style>
