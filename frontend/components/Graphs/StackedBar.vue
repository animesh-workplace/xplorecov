<template>
	<div>
		<Skeleton class="mb-2" height="24rem" v-if="isLoading" />
		<VChart ref="chart" class="w-[22rem] h-96" :option="graph_options" :class="{ hidden: isLoading }" />
	</div>
</template>

<script setup>
const props = defineProps({
	rawData: { type: Object, required: true, default: {} },
})

const chart = templateRef('chart')
const categories = Object.keys(props.rawData?.data)
const types = Object.keys(props.rawData?.data[categories[0]])
const seriesData = types.map((type) => ({
	name: type,
	type: 'bar',
	stack: 'total',
	data: categories.map((state) => props.rawData?.data[state][type]),
}))

const graph_options = ref({
	animation: true,
	title: {
		text: props.rawData?.name,
		textStyle: {
			width: 350,
			lineHeight: 20,
			overflow: 'break',
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
		data: categories,
	},
	yAxis: { type: 'value' },
	series: seriesData,
})

const isLoading = ref(true)

onMounted(() => {
	nextTick(() => {
		setTimeout(() => {
			isLoading.value = false
		}, 2000)
	})
})
</script>

<style scoped></style>
