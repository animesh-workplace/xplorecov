/** @type {import('tailwindcss').Config} */
// const defaultTheme = require('tailwindcss/defaultTheme')
// const fontFamily = defaultTheme.fontFamily
// fontFamily['sans'] = ['Lexend Deca', 'system-ui']

module.exports = {
	theme: {},
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