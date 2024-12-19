<template>
	<div>
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

const themeStore = useThemeStore()
onMounted(() => {
	nextTick(() => {
		themeStore.loadTheme()
	})
})
</script>
