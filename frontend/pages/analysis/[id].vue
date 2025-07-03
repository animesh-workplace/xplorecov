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
					<Icon name="tabler:chevron-right" class="w-4 h-4 mt-2" />
				</li>
				<li>
					<span class="text-surface-900 dark:text-surface-0 leading-normal">{{ route.params.id }}</span>
				</li>
			</ul>

			<div class="flex items-start flex-col lg:justify-between lg:flex-row">
				<div>
					<div class="font-medium text-3xl text-surface-900 dark:text-surface-0">Analysis Results</div>
					<div class="flex gap-6 items-center text-surface-700 dark:text-surface-100 flex-wrap">
						<div class="flex align-center items-center mt-4">
							<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" class="w-5 h-5 mr-2">
								<g fill="none" stroke="currentColor" stroke-width="1.5">
									<path
										d="M5 8c0-2.828 0-4.243.879-5.121C6.757 2 8.172 2 11 2h2c2.828 0 4.243 0 5.121.879C19 3.757 19 5.172 19 8v8c0 2.828 0 4.243-.879 5.121C17.243 22 15.828 22 13 22h-2c-2.828 0-4.243 0-5.121-.879C5 20.243 5 18.828 5 16zm0-3.924c-.975.096-1.631.313-2.121.803C2 5.757 2 7.172 2 10v4c0 2.828 0 4.243.879 5.121c.49.49 1.146.707 2.121.803M19 4.076c.975.096 1.631.313 2.121.803C22 5.757 22 7.172 22 10v4c0 2.828 0 4.243-.879 5.121c-.49.49-1.146.707-2.121.803"
									/>
									<path stroke-linecap="round" d="M9 13h6M9 9h6m-6 8h3" />
								</g>
							</svg>

							<span>{{ my_analysis?.total_sequences }} sequences</span>
						</div>

						<div class="flex align-center items-center mt-4">
							<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" class="w-7 h-7 mr-2">
								<path
									fill="currentColor"
									fill-rule="evenodd"
									d="M1.25 7A.75.75 0 0 1 2 6.25h8a.75.75 0 0 1 0 1.5H2A.75.75 0 0 1 1.25 7M17 7.75a4.25 4.25 0 1 0 0 8.5a4.25 4.25 0 0 0 0-8.5M11.25 12a5.75 5.75 0 1 1 11.5 0a5.75 5.75 0 0 1-11.5 0M17 9.25a.75.75 0 0 1 .75.75v1.566l.817.943a.75.75 0 0 1-1.134.982l-1-1.154a.75.75 0 0 1-.183-.49V10a.75.75 0 0 1 .75-.75M1.25 12a.75.75 0 0 1 .75-.75h6a.75.75 0 0 1 0 1.5H2a.75.75 0 0 1-.75-.75m0 5a.75.75 0 0 1 .75-.75h8a.75.75 0 0 1 0 1.5H2a.75.75 0 0 1-.75-.75"
									clip-rule="evenodd"
								/>
							</svg>

							<span>{{ analysis_duration }} mins duration</span>
						</div>

						<div class="flex align-center items-center mt-4">
							<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" class="w-6 h-6 mr-2">
								<g fill="none" stroke="currentColor" stroke-linecap="round" stroke-width="1.5">
									<path
										stroke-dasharray=".5 3.5"
										d="M22 12c0 5.523-4.477 10-10 10S2 17.523 2 12S6.477 2 12 2s10 4.477 10 10Z"
										opacity="0.5"
									/>
									<path d="M22 12c0-5.523-4.477-10-10-10" />
									<path stroke-linejoin="round" d="M12 9v4h4" />
								</g>
							</svg>

							<span>Expires {{ expires_in }}</span>
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

		<div class="grid grid-cols-6 gap-8 px-6 py-8 md:px-12 lg:px-20">
			<div class="col-span-2">
				<div class="w-full bg-white shadow-lg rounded-xl dark:bg-neutral-500">
					<div class="py-2 px-4 font-bold text-black text-md dark:text-white flex align-center">
						<svg
							viewBox="0 0 24 24"
							xmlns="http://www.w3.org/2000/svg"
							class="h-6 w-6 mr-2 stroke-slate-500 dark:stroke-slate-400 fill-none"
						>
							<path
								stroke-width="1.5"
								color="currentColor"
								stroke-linecap="round"
								stroke-linejoin="round"
								d="M3 4c0-1.655.345-2 2-2h4c1.655 0 2 .345 2 2s-.345 2-2 2H5c-1.655 0-2-.345-2-2m10 9c0-1.655.345-2 2-2h4c1.655 0 2 .345 2 2s-.345 2-2 2h-4c-1.655 0-2-.345-2-2m-9 7c0-1.655.345-2 2-2h4c1.655 0 2 .345 2 2s-.345 2-2 2H6c-1.655 0-2-.345-2-2m13-9c0-.465 0-.697-.038-.89a2 2 0 0 0-1.572-1.572c-.193-.038-.425-.038-.89-.038h-5c-.465 0-.697 0-.89-.038A2 2 0 0 1 7.038 6.89C7 6.697 7 6.465 7 6m10 9v1c0 1.886 0 2.828-.586 3.414S14.886 20 13 20h-1"
							/>
						</svg>

						<span> Workflow Steps </span>
					</div>
					<ul>
						<li
							:key="index"
							v-for="(step, index) in analysis_steps"
							:class="{ 'border-b-2': step.index != 7, 'border-t-2': step.index == 1 }"
							class="flex items-center justify-between py-2 text-gray-600 border-gray-100 dark:text-gray-100 dark:border-gray-800"
						>
							<div class="flex items-center justify-start text-sm">
								<span class="mx-4"> {{ step.index }} </span>
								<span class="mr-2"> {{ step.name }} </span>
								<span
									v-if="step.status == 'completed' && step.duration"
									class="text-xs text-gray-300"
								>
									(Took {{ step.duration }})
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

			<div class="col-span-4" v-if="analysis_complete">
				<div class="font-bold text-black text-md dark:text-white mt-2">Data Summary</div>
				<p v-for="(analysis, index) in my_analysis?.summary_reports" :key="index">
					{{ analysis.text_summary }}
				</p>
			</div>
		</div>

		<div class="mb-8 px-6 md:px-12 lg:px-20" v-if="analysis_complete">
			<Accordion value="">
				<AccordionPanel value="0" :pt="{ root: '!border-0' }">
					<AccordionHeader :pt="{ root: '!rounded-xl group/accordion' }">
						<div class="flex items-center gap-2">
							<Icon
								name="tabler:file-report"
								class="w-5 h-5 text-gray-400 group-hover/accordion:text-white"
							/>
							Combined analysis report
						</div>
					</AccordionHeader>
					<AccordionContent :pt="{ content: '!pt-4 !bg-[#121212]' }">
						<MyDataTable />
					</AccordionContent>
				</AccordionPanel>
			</Accordion>
		</div>

		<!-- Graphs -->
		<div class="grid grid-cols-4 gap-4 mb-8 px-6 md:px-12 lg:px-20" v-if="analysis_complete">
			<div v-for="(analysis, index) in my_analysis?.graph_reports" :key="index">
				<GraphsBar :rawData="analysis" v-if="analysis?.graph_type == 'Bar'" />
				<GraphsStackedBar :rawData="analysis" v-if="analysis?.graph_type == 'Stacked Bar'" />
			</div>
		</div>

		<!-- Chat Messages Container with proper scrolling -->

		<div
			class="space-y-3 text-center mt-24 mb-32"
			v-if="analysis_complete && my_analysis?.chat_messages?.length == 0"
		>
			<h3 class="text-gray-200 text-4xl font-semibold sm:text-5xl">Craving Deeper Insights?</h3>
			<p class="text-gray-400">
				Dive into your data with <strong>XPLORECoV-AI</strong> — ask complex questions, uncover hidden
				patterns, and instantly generate insightful visualizations.
			</p>
		</div>

		<div
			ref="chatContainer"
			class="h-96 mb-8 px-6 md:px-12 lg:px-20"
			v-if="analysis_complete && my_analysis?.chat_messages?.length > 0"
		>
			<div class="mt-24">
				<h3 class="text-gray-200 text-4xl font-semibold sm:text-5xl">Welcome to XPLORECoV-AI!</h3>
				<p class="text-gray-400">
					Unleash the power of AI to explore SARS-CoV-2 data. Ask anything about the virus, mutations, or
					sequences, and get instant insights and visualizations. Let’s get started!
				</p>
			</div>
			<div
				:key="index"
				:id="message.uuid"
				class="my-4 flex"
				v-for="(message, index) in my_analysis?.chat_messages"
				:class="{
					'justify-end': message.sender === 'human',
					'justify-start': message.sender === 'assistant',
					'pb-[calc(100vh-10rem)]': index === my_analysis.chat_messages.length - 1, // Only last message gets padding
				}"
			>
				<div
					class="max-w-full p-4 rounded-lg shadow-sm"
					:class="{
						'bg-blue-600 text-white rounded-tl-2xl rounded-tr-sm rounded-br-2xl rounded-bl-2xl':
							message.sender === 'human',
						'bg-gray-200 dark:bg-gray-700 text-gray-900 dark:text-white rounded-tr-2xl rounded-tl-sm rounded-bl-2xl rounded-br-2xl':
							message.sender === 'assistant',
					}"
				>
					<div v-if="message?.content == 'Loading'">
						<ProgressSpinner :pt="{ root: '!h-10 !w-10' }" stroke-width="4" />
					</div>
					<div v-else>
						<pre
							class="font-[Averta] text-wrap"
							v-if="message.sender === 'assistant'"
							v-typing="{
								typeSpeed: 10,
								hasCaret: false,
								text: message?.content,
							}"
						/>
						<div v-else v-html="message?.content" class="leading-relaxed" />
						<button
							type="button"
							v-if="message?.content.startsWith('Here are the two West Bengal genomes')"
							class="px-3 py-2 text-xs font-medium text-center text-white bg-blue-700 rounded-lg hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 dark:bg-blue-600 dark:hover:bg-blue-700 dark:focus:ring-blue-800"
						>
							Phylogenetic Tree
						</button>
					</div>
				</div>
			</div>
		</div>

		<div class="px-6 md:px-12 lg:px-20 fixed w-full bottom-5 shadow-2xl" v-if="analysis_complete">
			<div class="relative">
				<div class="absolute inset-y-0 start-0 flex items-center ps-3 pointer-events-none">
					<svg
						fill="none"
						aria-hidden="true"
						viewBox="0 0 20 20"
						xmlns="http://www.w3.org/2000/svg"
						class="w-4 h-4 text-gray-500 dark:text-gray-400"
					>
						<path
							stroke-width="2"
							stroke="currentColor"
							stroke-linecap="round"
							stroke-linejoin="round"
							d="m19 19-4-4m0-7A7 7 0 1 1 1 8a7 7 0 0 1 14 0Z"
						/>
					</svg>
				</div>
				<input
					type="search"
					v-model="search_prompt"
					@keypress.enter="AISearchQuery"
					class="block w-full p-4 pr-24 ps-10 text-sm text-gray-900 border border-gray-300 rounded-lg bg-gray-50 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white"
					placeholder="Ask your question to XPLORECoV-AI"
				/>
				<button
					disabled
					v-if="search_loading"
					class="absolute end-0 bottom-2 py-2 px-2.5 me-2 text-sm font-medium text-gray-900 bg-white rounded-lg border border-gray-200 hover:bg-gray-100 hover:text-blue-700 focus:z-10 focus:ring-4 focus:outline-none focus:ring-blue-700 focus:text-blue-700 dark:bg-gray-800 dark:text-gray-400 dark:border-gray-600 dark:hover:text-white dark:hover:bg-gray-700 inline-flex items-center"
				>
					<svg
						fill="none"
						role="status"
						aria-hidden="true"
						viewBox="0 0 100 101"
						xmlns="http://www.w3.org/2000/svg"
						class="inline w-4 h-4 mr-2 text-gray-200 animate-spin dark:text-gray-600"
					>
						<path
							fill="currentColor"
							d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z"
						/>
						<path
							fill="#1C64F2"
							d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z"
						/>
					</svg>
					Loading...
				</button>

				<button
					v-else
					@click="AISearchQuery"
					class="text-white absolute end-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 font-medium rounded-lg text-sm px-4 py-2 dark:bg-blue-600 dark:hover:bg-blue-700"
				>
					Search
				</button>
			</div>
		</div>
	</section>
</template>

<script setup>
import { round } from 'lodash-es'
import { v4 as uuidv4 } from 'uuid'
import { useUserAnalysis } from '@/api/analysis'
import { useSessionStore } from '@/stores/session'

const { data: my_analysis, error } = useAsyncData('specific_analysis', async () => {
	const route = useRoute()
	const { getSpecificAnalysis } = useUserAnalysis()
	return await getSpecificAnalysis(route.params.id)
})
const chatContainer = ref(null)
const dayjs = useDayjs()
const search_prompt = ref('')
const search_loading = ref(false)
const analysis_duration = computed(() => {
	const completionDate = my_analysis.value?.completion_date
	const submissionDate = my_analysis.value?.submission_date
	if (!completionDate || !submissionDate) return 0 // Return 0 if either is missing
	return round(dayjs(completionDate).diff(dayjs(submissionDate), 'minutes', true), 2)
})
const expires_in = computed(() => {
	const expirationDate = my_analysis.value?.expiration_date
	return dayjs().to(dayjs(expirationDate))
})

const scrollToId = (id) => {
	if (import.meta.client) {
		const element = document.getElementById(id)
		if (element) {
			const top = element.offsetTop
			window.scrollTo({ top, behavior: 'smooth' })
		}
	}
}

const AISearchQuery = async () => {
	if (!search_loading.value) {
		if (!search_prompt.value.trim()) return
		search_loading.value = true
		try {
			const saved_search_question = search_prompt.value
			// Add user message
			if (!my_analysis.value.chat_messages) {
				my_analysis.value.chat_messages = []
			}

			const scroll_id = uuidv4()
			my_analysis.value.chat_messages.push({
				content: search_prompt.value,
				content_type: 'text',
				created_at: new Date().toISOString(),
				parent_message_uuid: '',
				sender: 'human',
				uuid: scroll_id,
			})

			await new Promise((resolve) => setTimeout(resolve, 200))
			scrollToId(scroll_id)
			search_prompt.value = ''

			// Add assistant loading message
			my_analysis.value.chat_messages.push({
				content: 'Loading',
				content_type: 'text',
				created_at: new Date().toISOString(),
				parent_message_uuid: '',
				sender: 'assistant',
				uuid: uuidv4(),
			})

			// Simulate API call
			await new Promise((resolve) => setTimeout(resolve, 5000))

			// Replace loading message with actual response
			let content = ''
			if (
				saved_search_question == 'What are the most frequent spike protein mutations across all samples?'
			) {
				content = `Based on the combined analysis of 185 SARS-CoV-2 genome samples, the following are the top 10 most frequent spike (S) protein mutations observed:

				| Rank | Mutation     | Count | Frequency (%) |
				| ---- | ------------ | ----- | ------------- |
				| 1    | **S\:D614G** | 167   | **90.3%**     |
				| 2    | **S\:H655Y** | 163   | 88.1%         |
				| 3    | **S\:D796Y** | 160   | 86.5%         |
				| 4    | **S\:N679K** | 159   | 85.9%         |
				| 5    | **S\:Q954H** | 158   | 85.4%         |
				| 6    | **S\:N764K** | 156   | 84.3%         |
				| 7    | **S\:N969K** | 155   | 83.8%         |
				| 8    | **S\:G142D** | 149   | 80.5%         |
				| 9    | **S\:K417N** | 130   | 70.3%         |
				| 10   | **S\:N440K** | 124   | 67.0%         |
				`
			} else if (
				saved_search_question == 'What is the average number of substitutions and deletions per sample?'
			) {
				content = `Based on the combined analysis of the 185 SARS-CoV-2 genome samples:
				---

				🔹Average Number of Substitutions per Sample:
					Mean: 61.7 substitutions

				🔹 Average Number of Deletions per Sample:
					Mean: 3.2 deletions

				---

				Interpretation:

				* On average, each SARS-CoV-2 genome differs from the reference genome by ~62 point mutations (substitutions).
				* Each genome also carries about 3 deletions, which may impact protein function depending on their location and context.

				These values provide insight into the mutational burden across the viral population and are consistent with Omicron-lineage genomes, which are known to accumulate high numbers of substitutions, especially in the spike protein.`
			} else if (
				saved_search_question ==
				'Identify all SARS-CoV-2 sequences from West Bengal that contain both E484K and L452R mutations, and analyze their phylogenetic relationships to determine if they belong to emerging recombinant lineages distinct from Omicron and Delta.'
			) {
				content = `The result shows that zero West Bengal genomes in this report carry both mutations simultaneously. A quick scan of aaSubstitutions revealed many records with E484K plus other L452 changes (e.g. S:L452W, S:L452Q) but none with the exact L452R substitution.

Implications

Because no West Bengal sequences combine E484K + L452R in this dataset, there is no candidate set to build a phylogeny for a potential recombinant distinct from Omicron or Delta.
Any lineages displaying E484K with a different L452 mutation likely belong to known Omicron sub-lineages (e.g. XBB derivatives frequently show E484K with L452W) rather than an L452R-bearing Delta recombinant.`
			} else if (saved_search_question == 'do the same for L452W') {
				content = `Here are the two West Bengal genomes carrying the paired spike mutations $$S:E484K$$ and $$S:L452W$$:

{outputs_dict['seq_table']}

Both samples—INCOV293907 and INCOV293910a—were collected within the same week and sit in an Omicron-derived clade (XBB family).  Their shared signature (E484K + L452W) is unusual among West Bengal submissions, suggesting a common ancestor or a very recent local transmission chain rather than separate importations. They share the characteristic Spike mutations $$S:E484K$$ and $$S:L452W$$—a rare combination in West Bengal—alongside the typical XBB “backbone” (e.g.\ $$S:G142D$$, $$ORF9b:P10S$$, $$ORF1b:P314L$$).

Key points  
- Lineage calls: both are assigned to XBB.1.16-like sub-lineage bearing characteristic ORF1a:P3593F, ORF9b:P10S, etc.  
- Nextclade places them on the same terminal branch, 0–1 SNPs apart, reinforcing the notion of epidemiological linkage.  
- No Delta markers (e.g., $$ORF1b:P314L$$ without accompanying Delta spike set) appear, so they are not Delta–Omicron recombinants; instead they are part of the ongoing diversification of XBB.
`
			} else {
				content =
					'Here is the search result based on your query. This is a sample response that demonstrates the chat functionality.'
			}
			const lastMessageIndex = my_analysis.value.chat_messages.length - 1
			my_analysis.value.chat_messages[lastMessageIndex] = {
				...my_analysis.value.chat_messages[lastMessageIndex],
				content: content,
			}

			// Uncomment and modify this section for real API integration
			// try {
			// 	const { session } = useSessionStore()
			// 	const { askXPLORECoVAI } = useUserAnalysis()
			// 	await askXPLORECoVAI({
			// 		user_id: session,
			// 		analysis_id: route.params.id,
			// 		message: {
			// 			sender: 'human',
			// 			content_type: 'text',
			// 			parent_message_uuid: null,
			// 			content: search_prompt.value,
			// 		},
			// 	})
			// 	refreshAll()
			// } catch (err) {
			// 	console.log(err)
			// }
		} catch (err) {
			console.log(err)
		} finally {
			search_loading.value = false
		}
	}
}

// const AISearchQuery = async () => {
// 	search_loading.value = true
// 	try {
// 		my_analysis.value.chat_messages.push({
// 			content: search_prompt.value,
// 			content_type: '',
// 			created_at: '',
// 			parent_message_uuid: '',
// 			sender: 'human',
// 			uuid: '',
// 		})
// 		my_analysis.value.chat_messages.push({
// 			content: 'searching for result',
// 			content_type: '',
// 			created_at: '',
// 			parent_message_uuid: '',
// 			sender: 'assistant',
// 			uuid: '',
// 		})
// 		await new Promise((resolve) => setTimeout(resolve, 2000))
// 		scrollToBottom()
// 		// try {
// 		// 	const { session } = useSessionStore()
// 		// 	const { askXPLORECoVAI } = useUserAnalysis()
// 		// 	await askXPLORECoVAI({
// 		// 		user_id: session,
// 		// 		analysis_id: route.params.id,
// 		// 		message: {
// 		// 			sender: 'human',
// 		// 			content_type: 'text',
// 		// 			parent_message_uuid: null,
// 		// 			content: search_prompt.value,
// 		// 		},
// 		// 	})
// 		// 	refreshAll()
// 		// } catch (err) {
// 		// 	console.log(err)
// 		// }
// 	} catch (err) {
// 		console.log(err)
// 	}
// 	search_prompt.value = ''
// 	search_loading.value = false
// }

const route = useRoute()
const wsUrl = `ws://10.10.6.80/xplorecov/ws/analysis/${useCookie('session').value}/${route.params.id}/`
const { analysis_steps, tools_version, analysis_complete, disconnect } = useWebSocket(wsUrl)

const refreshAll = async () => {
	try {
		const { getSpecificAnalysis } = useUserAnalysis()
		my_analysis.value = await getSpecificAnalysis(route.params.id)
		disconnect(my_analysis.value?.overall_status == 'SUCCESS')
	} catch (err) {
		console.log(err)
	}
}

watch(
	() => analysis_complete.value,
	async (newValue) => {
		if (newValue) {
			refreshAll()
		}
	},
	{ immediate: true },
)

const wsChatUrl = `ws://10.10.6.80/xplorecov/ws/chat/${useCookie('session').value}/${route.params.id}/`
const { message, Chatdisconnect } = useChatWebSocket(wsChatUrl)
const { displayedText } = useTyping(message.value?.message, 30)

watch(
	() => message.value?.message,
	async (newValue) => {
		if (newValue) {
			refreshAll()
		}
	},
	{ immediate: true },
)
</script>

<style scoped>
.cursor {
	animation: blink 0.8s infinite;
}

@keyframes blink {
	50% {
		opacity: 0;
	}
}

.spinner {
	animation: spin 1s linear infinite;
}

@keyframes spin {
	from {
		transform: rotate(0deg);
	}
	to {
		transform: rotate(360deg);
	}
}
</style>
