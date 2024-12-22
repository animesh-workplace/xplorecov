import withNuxt from './.nuxt/eslint.config.mjs'

export default withNuxt({
	files: ['**/*.vue', '***/*.js'],
	rules: {
		'vue/attributes-order': 'off',
		'vue/no-v-model-argument': 'off',
		semi: ['error', 'never'],
	},
})
