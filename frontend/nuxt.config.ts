import Aura from '@primevue/themes/aura'

export default defineNuxtConfig({
	app: { head: { title: 'XPLORECoV | National Institute of Biomedical Genomics', meta: [], link: [] } },
	modules: [
		'dayjs-nuxt',
		'nuxt-umami',
		'@nuxt/icon',
		'@pinia/nuxt',
		'@nuxt/fonts',
		'@nuxt/image',
		'nuxt-echarts',
		'@nuxt/eslint',
		'@vueuse/nuxt',
		'notivue/nuxt',
		'@nuxtjs/device',
		'@nuxtjs/tailwindcss',
		'@primevue/nuxt-module',
	],
	echarts: {
		ssr: true,
		charts: ['BarChart'],
		components: ['GridComponent', 'DatasetComponent', 'TooltipComponent', 'ToolboxComponent'],
	},
	devtools: { enabled: true },
	compatibilityDate: '2024-11-01',
	notivue: { position: 'bottom-right' },
	primevue: { options: { ripple: true, theme: { preset: Aura, options: { darkModeSelector: '.app-dark' } } } },
	css: [
		'@/assets/css/main.css',
		'primeicons/primeicons.css',
		'notivue/notification.css',
		'notivue/animations.css',
	],
	umami: { enabled: false, id: '12d666c0-d0bf-4271-ae3b-2ff52f81be58', host: 'https://research.nibmg.ac.in' },
	dayjs: { plugins: ['customParseFormat'] },
})
