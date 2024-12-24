import withNuxt from './.nuxt/eslint.config.mjs'

export default withNuxt({
	files: ['**/*.vue', '***/*.js'],
	rules: {
		'vue/attributes-order': 'off',
		'vue/valid-v-model': 'off',
		semi: ['error', 'never'],
	},
})
