import Aura from '@primevue/themes/aura'
import { definePreset } from '@primevue/themes'

const Xaura = definePreset(Aura, {
	semantic: {
		colorScheme: {
			dark: {
				surface: {
					0: '#ffffff',
					50: '{neutral.50}',
					100: '{neutral.100}',
					200: '{neutral.200}',
					300: '{neutral.300}',
					400: '{neutral.400}',
					500: '{neutral.500}',
					600: '{neutral.600}',
					700: '{neutral.700}',
					800: '{neutral.800}',
					900: '{neutral.700}',
					950: '{neutral.950}',
				},
			},
		},
	},
})

export default defineNuxtConfig({
<<<<<<< HEAD
=======
	runtimeConfig: {
		public: {
			API_BASE_URL: process.env.API_BASE_URL || 'http://localhost:8009/api',
		},
	},
>>>>>>> main#1
	app: {
		head: {
			title: 'XPLORECoV | National Institute of Biomedical Genomics',
			meta: [],
			link: [],
		},
	},
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
	dayjs: { plugins: ['customParseFormat'] },
	umami: { enabled: false, id: '12d666c0-d0bf-4271-ae3b-2ff52f81be58', host: 'https://research.nibmg.ac.in' },
	primevue: { options: { ripple: true, theme: { preset: Xaura, options: { darkModeSelector: '.dark' } } } },
	css: [
		'primeicons/primeicons.css',
		'notivue/notification.css',
		'notivue/animations.css',
		'@/assets/css/main.css',
		'@/assets/css/fonts.css',
		'@/assets/css/jquery.mCustomScrollbar.css',
	],
	// build: { transpile: ['lodash'] },
})
