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
</script>
