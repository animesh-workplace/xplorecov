import { defineStore } from 'pinia'
import { useStorage } from '@vueuse/core'

export const useThemeStore = defineStore('theme', {
	state: () => ({ isDarkMode: false }),

	actions: {
		toggleTheme() {
			this.isDarkMode = !this.isDarkMode
			this.applyTheme()
			if (import.meta.client) {
				useStorage('theme', this.isDarkMode ? 'dark' : 'light')
			}
		},
		loadTheme() {
			if (import.meta.client) {
				const savedTheme = useStorage('theme')
				if (savedTheme.value) {
					this.isDarkMode = savedTheme.value === 'dark'
				}
				this.applyTheme()
			}
		},
		applyTheme() {
			const themeClass = 'app-dark'
			if (this.isDarkMode) {
				document.documentElement.classList.add(themeClass)
			} else {
				document.documentElement.classList.remove(themeClass)
			}
		},
	},
})
