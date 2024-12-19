<template>
	<div>
		<!-- <Button @click="toggleDarkMode" rounded size="large" label="Large button" icon="pi pi-check" /> -->
		<Button
			rounded
			severity="secondary"
			@click="toggleDarkMode"
			:icon="defaultTheme ? 'pi pi-sun' : 'pi pi-moon'"
		/>
		<!-- <div v-if="isLoading" class="loader">Loading Graph...</div> -->
		<!-- <VChart ref="chart" class="test" :option="option" :class="{ hidden: isLoading }" /> -->

		<div class="mx-8 my-8 grid lg:grid-cols-2 lg:gap-14 grid-cols-1">
			<UploadMetadata @verification_status="verifyMetadata" />
			<UploadSequence @verification_status="verifySequence" />
		</div>
	</div>
</template>

<script setup lang="ts">
import { ref, onMounted, nextTick } from 'vue'
const chart = templateRef('chart')

const defaultTheme = ref(false)
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
	{ name: 'Missing Metadata/Sequence check', verification: false, data: [], error: null, show: false },
	{ name: 'Duplicate check', verification: false, data: [], error: null, show: false },
	{ name: 'Already present check', verification: false, data: [], error: null, show: false },
])

function toggleDarkMode() {
	defaultTheme.value = !defaultTheme.value
	document.documentElement.classList.toggle('app-dark')
}

const verifyMetadata = (metadata) => {
	all_qc_checks.value[0].verification = metadata.verification
	all_qc_checks.value[1].verification = metadata.verification
}

const verifySequence = (sequence) => {
	all_qc_checks.value[2].verification = sequence.verification
	all_qc_checks.value[3].verification = sequence.verification
}

// function random() {
// 	return Math.round(300 + Math.random() * 700) / 10
// }

// function getData(): ECOption {
// 	return {
// 		animation: true,
// 		tooltip: {
// 			className: 'echarts-tooltip',
// 		},
// 		toolbox: {
// 			show: true,
// 			feature: {
// 				dataZoom: {},
// 				saveAsImage: {},
// 			},
// 		},
// 		dataset: {
// 			dimensions: ['Product', '2015', '2016', '2017'],
// 			source: [
// 				{
// 					Product: 'Matcha Latte',
// 					2015: random(),
// 					2016: random(),
// 					2017: random(),
// 				},
// 				{
// 					Product: 'Milk Tea',
// 					2015: random(),
// 					2016: random(),
// 					2017: random(),
// 				},
// 				{
// 					Product: 'Cheese Cocoa',
// 					2015: random(),
// 					2016: random(),
// 					2017: random(),
// 				},
// 				{
// 					Product: 'Walnut Brownie',
// 					2015: random(),
// 					2016: random(),
// 					2017: random(),
// 				},
// 			],
// 		},
// 		xAxis: { type: 'category' },
// 		yAxis: {},
// 		itemStyle: { borderRadius: [5, 5, 0, 0] },
// 		series: [{ type: 'bar' }, { type: 'bar' }, { type: 'bar' }],
// 	}
// }

// const isLoading = ref(true)

// const option = shallowRef(getData())
// function refreshData() {
// 	option.value = getData()
// }

onMounted(() => {
	nextTick(() => {
		// isLoading.value = false
	})
})
</script>

<style>
.test {
	margin: 5px;
	width: 700px;
	height: 500px;
}
</style>
