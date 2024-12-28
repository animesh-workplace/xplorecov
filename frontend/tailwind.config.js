/** @type {import('tailwindcss').Config} */
const defaultTheme = require('tailwindcss/defaultTheme')
const fontFamily = defaultTheme.fontFamily
fontFamily['sans'] = ['Averta', 'system-ui']

module.exports = {
	theme: {},
	// darkMode: ['class', 'app-dark'],
	content: [
		'./pages/**/*.vue',
		'./layouts/**/*.vue',
		'./nuxt.config.{js,ts}',
		'./plugins/**/*.{js,ts}',
		'./components/**/*.{js,vue,ts}',
		'./node_modules/flowbite.{js,ts}',
	],
	plugins: [],
}
