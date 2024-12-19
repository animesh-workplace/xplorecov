<template>
	<div>
		<Notivue v-slot="item">
			<Notification :item="item" :theme="materialTheme" />
		</Notivue>
		<NuxtPage />
	</div>
</template>

<script setup>
import { materialTheme } from 'notivue'
</script>
