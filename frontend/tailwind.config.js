const defaultTheme = require('tailwindcss/defaultTheme')

module.exports = {
	content: [
		'./pages/**/*.vue',
		'./layouts/**/*.vue',
		'./nuxt.config.{js,ts}',
		'./plugins/**/*.{js,ts}',
		'./components/**/*.{js,vue,ts}',
		'./node_modules/flowbite/**/*.js',
	],
	darkMode: 'class',
	theme: {
		extend: {
			fontFamily: {
				sans: ['Averta', ...defaultTheme.fontFamily.sans],
			},
		},
	},
	plugins: [require('tailwind-scrollbar')],
}
