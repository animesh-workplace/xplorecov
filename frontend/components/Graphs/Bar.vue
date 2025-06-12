<template>
	<div>
		<Skeleton class="mb-2" height="24rem" v-if="isLoading" />
		<VChart
			autoresize
			ref="chart"
			:option="graph_options"
			class="w-[22rem] h-[22rem]"
			:class="{ hidden: isLoading }"
		/>
	</div>
</template>

<script setup>
const props = defineProps({
	rawData: { type: Object, required: true, default: {} },
})

const chart = templateRef('chart')
const values = Object.values(props.rawData?.data)
const categories = Object.keys(props.rawData?.data)

const graph_options = ref({
	yAxis: {
		type: 'value',
		axisLabel: {
			fontFamily: 'Averta',
			fontWeight: 500,
		},
	},
	tooltip: {
		trigger: 'axis',
		position: 'right',
		renderMode: 'richText',
		axisPointer: { type: 'shadow' },
		borderColor: '#fff',
		textStyle: {
			fontFamily: 'Averta',
			fontWeight: 500,
		},
		order: 'valueDesc',
	},
	animation: true,
	title: {
		text: props.rawData?.name,
		textStyle: {
			width: 400,
			lineHeight: 20,
			overflow: 'break',
		},
	},
	grid: {
		left: '0%',
		right: '0%',
		bottom: '0%',
		containLabel: true,
	},
	series: [{ type: 'bar', data: values }],
	itemStyle: { borderRadius: [5, 5, 0, 0] },
	xAxis: {
		type: 'category',
		data: categories,
		axisLabel: {
			fontFamily: 'Averta',
			fontWeight: 500,
		},
	},
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
