import withNuxt from './.nuxt/eslint.config.mjs'

export default withNuxt({
	files: ['**/*.vue', '***/*.js'],
	rules: {
		'vue/attributes-order': 'off',
		semi: ['error', 'never'],
	},
})
