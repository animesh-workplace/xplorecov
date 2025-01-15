<template>
	<section>
		<div class="bg-surface-0 dark:bg-surface-950 px-6 py-8 md:px-12 lg:px-20">
			<ul class="list-none p-0 m-0 flex items-center font-medium mb-4">
				<li>
					<NuxtLink
						to="/analysis"
						class="text-surface-500 dark:text-surface-300 no-underline leading-normal cursor-pointer hover:text-slate-400"
					>
						Analysis
					</NuxtLink>
				</li>
				<li class="px-2">
					<i class="pi pi-angle-right text-surface-500 dark:text-surface-300 leading-normal" />
				</li>
				<li>
					<span class="text-surface-900 dark:text-surface-0 leading-normal">{{ route.params.id }}</span>
				</li>
			</ul>
			<div class="flex items-start flex-col lg:justify-between lg:flex-row">
				<div>
					<div class="font-medium text-3xl text-surface-900 dark:text-surface-0">Analysis Results</div>
					<div class="flex items-center text-surface-700 dark:text-surface-100 flex-wrap">
						<div class="mr-8 flex items-center mt-4">
							<i class="pi pi-users mr-2" />
							<span>332 Active Users</span>
						</div>
						<div class="mr-8 flex items-center mt-4">
							<i class="pi pi-globe mr-2" />
							<span>9402 Sessions</span>
						</div>
						<div class="flex items-center mt-4">
							<i class="pi pi-clock mr-2" />
							<span>2.32m Avg. Duration</span>
						</div>
					</div>
				</div>
			</div>
		</div>

		<div class="flex gap-4 justify-center px-6 md:px-12 lg:px-20">
			<div
				:key="index"
				v-for="(item, index) in tools_version.tools"
				class="relative w-64 py-3 px-4 overflow-hidden bg-white dark:bg-neutral-500 shadow-lg rounded-xl"
			>
				<svg
					viewBox="0 0 32 32"
					xmlns="http://www.w3.org/2000/svg"
					class="absolute w-20 h-20 mb-4 -right-8 -bottom-5 fill-slate-500 dark:fill-slate-400"
				>
					<path
						d="M25 21c-.74 0-1.424.216-2.02.567l-2.03-2.031l-1.415 1.414l2.032 2.031A3.95 3.95 0 0 0 21 25c0 2.206 1.794 4 4 4c.356 0 .694-.061 1.023-.15l-2.437-2.436A2 2 0 0 1 23 25a2.002 2.002 0 0 1 3.414-1.414l2.437 2.437A4 4 0 0 0 29 25c0-2.206-1.795-4-4-4m-4.05-8.536L24.714 8.7c.391.187.824.3 1.286.3a3 3 0 1 0-3-3c0 .462.113.895.3 1.286l-3.764 3.764zM26 5c.551 0 1 .449 1 1s-.449 1-1 1s-1-.449-1-1s.449-1 1-1m-10 7a4 4 0 0 0-4 4c0 .74.215 1.425.566 2.02l-5.28 5.28A3 3 0 0 0 6 23a3 3 0 1 0 3 3c0-.462-.113-.895-.3-1.286l5.28-5.28A3.96 3.96 0 0 0 16 20a4 4 0 0 0 0-8M6 27a1.001 1.001 0 0 1 0-2a1.001 1.001 0 0 1 0 2m10-9c-1.103 0-2-.897-2-2s.897-2 2-2s2 .897 2 2s-.897 2-2 2m-9-7c.74 0 1.424-.215 2.02-.567l2.03 2.031l1.414-1.414l-2.03-2.031C10.783 8.424 11 7.739 11 7c0-2.206-1.794-4-4-4c-.356 0-.694.062-1.023.15l2.437 2.436A2.002 2.002 0 0 1 7 9a2 2 0 0 1-1.414-.586L3.149 5.977A4 4 0 0 0 3 7c0 2.206 1.794 4 4 4"
					/>
				</svg>

				<div class="w-4/6">
					<p class="text-lg font-medium text-gray-800 dark:text-gray-100 capitalize">
						{{ index.replace('_version', '') }}
					</p>
					<p class="text-gray-500 dark:text-gray-200">{{ item }}</p>
				</div>
			</div>
		</div>

		<div class="px-6 py-8 md:px-12 lg:px-20 lg:max-w-screen-sm">
			<div class="w-full bg-white shadow-lg rounded-xl dark:bg-neutral-500">
				<p class="py-2 px-4 font-bold text-black text-md dark:text-white">Analysis Steps</p>
				<ul>
					<li
						:key="index"
						v-for="(step, index) in analysis_steps"
						:class="{ 'border-b-2': step.index != 6, 'border-t-2': step.index == 1 }"
						class="flex items-center justify-between py-2 text-gray-600 border-gray-100 dark:text-gray-100 dark:border-gray-800"
					>
						<div class="flex items-center justify-start text-sm">
							<span class="mx-4"> {{ step.index }} </span>
							<span class="mr-2"> {{ step.name }} </span>
							<span v-if="step.status == 'completed' && step.duration" class="text-xs text-gray-300">
								(Took {{ step.duration }} secs)
							</span>
						</div>

						<i
							v-if="step.status == 'pending'"
							class="pi pi-pause-circle mx-4 text-gray-400 dark:text-gray-300"
						/>

						<i
							v-if="step.status == 'completed'"
							class="pi pi-check-circle mx-4 text-green-400 dark:text-green-500"
						/>

						<svg
							width="17"
							height="17"
							viewBox="0 0 24 24"
							v-if="step.status == 'loading'"
							xmlns="http://www.w3.org/2000/svg"
							class="mx-4 text-gray-400 dark:text-gray-300 stroke-yellow-700"
						>
							<g class="spinner">
								<circle cx="12" cy="12" r="9.5" fill="none" stroke-width="3" />
							</g>
						</svg>
					</li>
				</ul>
			</div>
		</div>
	</section>
</template>

<script setup>
// import { useUserAnalysis } from '@/api/analysis'

// const { data: my_analysis, error } = useAsyncData('analysis', async () => {
// 	const { getAnalysis } = useUserAnalysis()
// 	return await getAnalysis()
// })

const route = useRoute()
const wsUrl = `ws://10.10.6.80/xplorecov/ws/analysis/${useCookie('session').value}/${route.params.id}/`
const { analysis_steps, tools_version } = useWebSocket(wsUrl)
</script>

<style scoped></style>
