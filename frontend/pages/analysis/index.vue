<template>
	<div>
		<div class="bg-surface-0 dark:bg-surface-950 px-6 py-8 md:px-12 lg:px-20">
			<ul class="list-none p-0 m-0 flex items-center font-medium mb-4">
				<li>
					<a class="text-surface-500 dark:text-surface-300 no-underline leading-normal cursor-pointer">
						Application
					</a>
				</li>
				<li class="px-2">
					<i class="pi pi-angle-right text-surface-500 dark:text-surface-300 leading-normal" />
				</li>
				<li>
					<span class="text-surface-900 dark:text-surface-0 leading-normal">Analytics</span>
				</li>
			</ul>
		</div>

		<div class="px-6 py-8 md:px-12 lg:px-20">
			<DataTable
				rowHover
				paginator
				:rows="10"
				size="small"
				:value="my_analysis"
				:rowsPerPageOptions="[5, 10, 20, 50]"
				:pt="{ pcPaginator: { root: '!rounded-b-lg' } }"
			>
				<Column
					field="index"
					header="Sl. No"
					class="!text-center"
					:pt="{ columnHeaderContent: '!justify-center', headerCell: '!rounded-tl-lg !p-4' }"
				/>
				<Column
					header="Analysis ID"
					class="!text-center"
					field="analysis_id"
					:pt="{ columnHeaderContent: '!justify-center' }"
				/>
				<Column
					class="!text-center"
					field="submission_date"
					header="Submission Date"
					:pt="{ columnHeaderContent: '!justify-center' }"
				/>
				<Column
					class="!text-center"
					field="total_sequences"
					header="Total Sequences"
					:pt="{ columnHeaderContent: '!justify-center' }"
				/>
				<Column
					header="Status"
					class="!text-center"
					field="overall_status"
					:pt="{ columnHeaderContent: '!justify-center' }"
				>
					<template #body="{ data }">
						<Tag
							:value="data.overall_status"
							:severity="
								data.overall_status == 'SUCCESS'
									? 'success'
									: data.overall_status == 'ERROR'
										? 'danger'
										: 'warn'
							"
						/>
					</template>
				</Column>
				<Column class="w-44 !text-center" :pt="{ headerCell: '!rounded-tr-lg' }">
					<template #body="{ data }">
						<NuxtLink :to="`/analysis/${data.analysis_id}`">
							<Button
								size="small"
								label="View Analysis"
								severity="secondary"
								icon="pi pi-external-link"
							/>
						</NuxtLink>
					</template>
				</Column>
			</DataTable>
		</div>
	</div>
</template>

<script setup>
import { map } from 'lodash'
import { useUserAnalysis } from '@/api/analysis'

const { data: my_analysis, error } = useAsyncData('analysis', async () => {
	const dayjs = useDayjs()
	const { getAnalysis } = useUserAnalysis()
	const response = await getAnalysis()

	return map(response, (d, i) => ({
		...d,
		index: i + 1,
		submission_date: dayjs(d.submission_date.replace(/\.\d{3}\d*/, '')).format('DD-MM-YYYY hh:mm a'),
	}))
})
</script>

<style scoped></style>
