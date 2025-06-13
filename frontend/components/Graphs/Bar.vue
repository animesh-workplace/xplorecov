<template>
	<div>
		<Skeleton class="mb-2" height="24rem" v-if="isLoading" />
		<div v-show="!isLoading" class="h-[25rem] w-full">
			<VChart ref="chart" :option="graph_options" class="w-full h-full" autoresize />
		</div>
	</div>
</template>

<script setup>
const props = defineProps({
	rawData: { type: Object, required: true, default: () => ({}) },
})

const isLoading = ref(true)
const chart = templateRef('chart')

const graph_options = ref({
	animation: true,
	yAxis: {
		type: 'value',
		axisLabel: {
			fontFamily: 'Averta',
			fontWeight: 500,
		},
	},
	xAxis: {
		type: 'category',
		data: [],
		axisLabel: {
			fontFamily: 'Averta',
			fontWeight: 500,
		},
	},
	tooltip: {
		trigger: 'axis',
		position: 'right',
		axisPointer: { type: 'shadow' },
		borderColor: '#fff',
		textStyle: {
			fontFamily: 'Averta',
			fontWeight: 500,
		},
		order: 'valueDesc',
	},
	title: {
		text: '',
		textStyle: {
			width: 400,
			lineHeight: 20,
			overflow: 'break',
		},
	},
	grid: {
		top: '15%', // Give space for title
		left: '0%',
		right: '4%',
		bottom: '0%',
		containLabel: true,
	},
	series: [
		{
			type: 'bar',
			data: [],
			itemStyle: {
				borderRadius: [5, 5, 0, 0],
			},
		},
	],
})

const updateChart = () => {
	if (!props.rawData?.data) {
		console.warn('No data available to update chart')
		return
	}

	// Update the options
	graph_options.value = {
		...graph_options.value,
		title: {
			...graph_options.value.title,
			text: props.rawData?.name || '',
		},
		xAxis: {
			...graph_options.value.xAxis,
			data: Object.keys(props.rawData.data),
		},
		series: [
			{
				...graph_options.value.series[0],
				data: Object.values(props.rawData.data),
			},
		],
	}

	// Force chart to re-render if it exists
	nextTick(() => {
		if (chart.value) {
			chart.value.setOption(graph_options.value, true)
		}
	})
}

// Watch for data changes
watch(
	() => props.rawData,
	(newValue) => {
		if (newValue?.data && Object.keys(newValue.data).length > 0) {
			updateChart()
			isLoading.value = false
		}
	},
	{ deep: true, immediate: true },
)

const initializeChart = () => {
	// Wait for DOM to be ready and have dimensions
	setTimeout(() => {
		if (chart.value && chart.value.$el) {
			const chartElement = chart.value.$el
			if (chartElement.clientWidth > 0 && chartElement.clientHeight > 0) {
				if (props.rawData?.data && Object.keys(props.rawData.data).length > 0) {
					updateChart()
				}
				isLoading.value = false
			} else {
				// Retry if dimensions are still 0
				console.warn('Chart dimensions still 0, retrying...')
				initializeChart()
			}
		}
	}, 100)
}

onMounted(() => {
	nextTick(() => {
		initializeChart()
	})
})
</script>
