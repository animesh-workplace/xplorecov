<template>
	<div>
		<Skeleton class="mb-2" height="24rem" v-if="isLoading" />
		<div v-show="!isLoading" class="h-[25rem] w-full">
			<VChart autoresize ref="chart" :option="graph_options" class="w-full h-full" />
		</div>
	</div>
</template>

<script setup>
const props = defineProps({
	rawData: { type: Object, required: true, default: () => ({}) },
})

const isLoading = ref(true)
const chart = templateRef('chart')

// Make data processing reactive
const categories = computed(() => {
	return props.rawData?.data ? Object.keys(props.rawData.data) : []
})

const types = computed(() => {
	if (!props.rawData?.data || categories.value.length === 0) return []
	return Object.keys(props.rawData.data[categories.value[0]] || {})
})

const seriesData = computed(() => {
	if (!props.rawData?.data || categories.value.length === 0 || types.value.length === 0) {
		return []
	}

	return types.value.map((type) => ({
		name: type,
		type: 'bar',
		stack: 'total',
		data: categories.value.map((state) => props.rawData.data[state]?.[type] || 0),
	}))
})

const graph_options = computed(() => ({
	animation: true,
	title: {
		text: props.rawData?.name || '',
		textStyle: {
			width: 350,
			lineHeight: 20,
			fontWeight: 500,
			overflow: 'break',
			fontFamily: 'Averta',
		},
	},
	grid: {
		top: '15%', // Give space for title
		left: '0%',
		right: '4%',
		bottom: '0%',
		containLabel: true,
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
		formatter: function (params) {
			let tooltipContent = `${params[0].axisValue}<br/>`
			params.forEach((item) => {
				if (item.data > 0) {
					tooltipContent += `${item.marker} ${item.seriesName}: ${item.data}<br/>`
				}
			})
			return tooltipContent
		},
	},
	xAxis: {
		type: 'category',
		data: categories.value,
	},
	yAxis: { type: 'value' },
	series: seriesData.value,
}))

// Watch for data changes and update loading state
watch(
	() => props.rawData,
	(newData) => {
		if (newData && newData.data && Object.keys(newData.data).length > 0) {
			isLoading.value = false
		}
	},
	{ immediate: true, deep: true },
)

onMounted(() => {
	nextTick(() => {
		// Check if data is already available
		if (props.rawData && props.rawData.data && Object.keys(props.rawData.data).length > 0) {
			isLoading.value = false
		}

		// Ensure chart resizes properly after loading state changes
		setTimeout(() => {
			if (chart.value && chart.value.resize) {
				chart.value.resize()
			}
		}, 100)
	})
})
</script>

<style scoped></style>
