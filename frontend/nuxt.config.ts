import Aura from '@primevue/themes/aura'

export default defineNuxtConfig({
	compatibilityDate: '2024-11-01',
	devtools: { enabled: true },
	modules: [
		'dayjs-nuxt',
		'nuxt-umami',
		'@nuxt/icon',
		'@pinia/nuxt',
		'@nuxt/fonts',
		'@nuxt/image',
		'nuxt-lodash',
		'nuxt-echarts',
		'@nuxt/eslint',
		'@vueuse/nuxt',
		'@nuxtjs/device',
		'@nuxtjs/tailwindcss',
		'@primevue/nuxt-module',
	],
	primevue: {
		options: {
			ripple: true,
			theme: { preset: Aura },
		},
	},
	umami: { enabled: false, id: '12d666c0-d0bf-4271-ae3b-2ff52f81be58', host: 'https://research.nibmg.ac.in' },
	echarts: {
		ssr: true,
		charts: ['BarChart'],
		components: ['GridComponent', 'DatasetComponent', 'TooltipComponent', 'ToolboxComponent'],
	},
})
