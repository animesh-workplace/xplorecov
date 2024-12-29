<template>
	<div class="min-h-screen">
		<Notivue v-slot="item">
			<Notification :item="item" :theme="materialTheme" />
		</Notivue>

		<NuxtLayout>
			<NuxtPage />
		</NuxtLayout>
	</div>
</template>

<script setup>
import { materialTheme } from 'notivue'
import { useThemeStore } from '@/stores/theme'
import { useSessionStore } from '@/stores/session'

const themeStore = useThemeStore()
const sessionStore = useSessionStore()

onMounted(() => {
	themeStore.loadTheme()
	sessionStore.getSessionInfo()
})
</script>
