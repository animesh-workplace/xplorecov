<template>
	<div class="card">
		<DataTable
			paginator
			:rows="10"
			dataKey="id"
			:value="customers"
			:loading="loading"
			v-model:filters="filters"
			:pt="{
				headerRow: '!py-4',
				header: 'rounded-xl mb-1',
				tableContainer: 'rounded-xl',
				column: { headerCell: '!py-4' },
				pcPaginator: { paginatorContainer: '!border-0 pt-1', root: '!rounded-xl' },
			}"
		>
			<template #empty> No data found. </template>
			<template #loading> Loading data. Please wait. </template>

			<Column field="name" header="Name">
				<template #body="{ data }">
					{{ data.name }}
				</template>
			</Column>

			<Column header="Country" filterField="country.name">
				<template #body="{ data }">
					<div class="flex items-center gap-2">
						<img
							alt="flag"
							src="https://primefaces.org/cdn/primevue/images/flag/flag_placeholder.png"
							:class="`flag flag-${data.country.code}`"
							style="width: 24px"
						/>
						<span>{{ data.country.name }}</span>
					</div>
				</template>
			</Column>

			<Column header="Agent" filterField="representative" :showFilterMenu="false">
				<template #body="{ data }">
					<div class="flex items-center gap-2">
						<img
							:alt="data.representative.name"
							:src="`https://primefaces.org/cdn/primevue/images/avatar/${data.representative.image}`"
							style="width: 32px"
						/>
						<span>{{ data.representative.name }}</span>
					</div>
				</template>
			</Column>

			<Column field="status" header="Status" :showFilterMenu="false">
				<template #body="{ data }">
					<Tag :value="data.status" :severity="getSeverity(data.status)" />
				</template>
			</Column>

			<Column field="verified" header="Verified" dataType="boolean">
				<!-- <template #body="{ data }">
					<i
						class="pi"
						:class="{
							'pi-check-circle text-green-500': data.verified,
							'pi-times-circle text-red-400': !data.verified,
						}"
					></i>
				</template> -->
			</Column>
		</DataTable>
	</div>
</template>

<script setup>
import { FilterMatchMode } from '@primevue/core/api'

const customers = ref()
const filters = ref({
	global: { value: null, matchMode: FilterMatchMode.CONTAINS },
	name: { value: null, matchMode: FilterMatchMode.STARTS_WITH },
	'country.name': { value: null, matchMode: FilterMatchMode.STARTS_WITH },
	representative: { value: null, matchMode: FilterMatchMode.IN },
	status: { value: null, matchMode: FilterMatchMode.EQUALS },
	verified: { value: null, matchMode: FilterMatchMode.EQUALS },
})
const representatives = ref([
	{ name: 'Amy Elsner', image: 'amyelsner.png' },
	{ name: 'Anna Fali', image: 'annafali.png' },
	{ name: 'Asiya Javayant', image: 'asiyajavayant.png' },
	{ name: 'Bernardo Dominic', image: 'bernardodominic.png' },
	{ name: 'Elwin Sharvill', image: 'elwinsharvill.png' },
	{ name: 'Ioni Bowcher', image: 'ionibowcher.png' },
	{ name: 'Ivan Magalhaes', image: 'ivanmagalhaes.png' },
	{ name: 'Onyama Limba', image: 'onyamalimba.png' },
	{ name: 'Stephen Shaw', image: 'stephenshaw.png' },
	{ name: 'XuXue Feng', image: 'xuxuefeng.png' },
])
const statuses = ref(['unqualified', 'qualified', 'new', 'negotiation', 'renewal', 'proposal'])
const loading = ref(true)

onMounted(() => {
	customers.value = [
		{
			id: 1000,
			name: 'James Butt',
			country: { name: 'Algeria', code: 'dz' },
			company: 'Benton, John B Jr',
			date: '2015-09-13',
			status: 'unqualified',
			verified: true,
			activity: 17,
			representative: { name: 'Ioni Bowcher', image: 'ionibowcher.png' },
			balance: 70663,
		},
		{
			id: 1001,
			name: 'Josephine Darakjy',
			country: { name: 'Argentina', code: 'ar' },
			company: 'Chanay, Jeffrey A Esq',
			date: '2019-01-05',
			status: 'qualified',
			verified: false,
			activity: 42,
			representative: { name: 'Amy Elsner', image: 'amyelsner.png' },
			balance: 45890,
		},
		{
			id: 1002,
			name: 'Art Venere',
			country: { name: 'Brazil', code: 'br' },
			company: 'Chemel, James L Cpa',
			date: '2020-07-21',
			status: 'new',
			verified: true,
			activity: 12,
			representative: { name: 'Asiya Javayant', image: 'asiyajavayant.png' },
			balance: 11890,
		},
		{
			id: 1003,
			name: 'Lenna Paprocki',
			country: { name: 'Canada', code: 'ca' },
			company: 'Feltz Printing Service',
			date: '2018-11-02',
			status: 'negotiation',
			verified: false,
			activity: 54,
			representative: { name: 'Onyama Limba', image: 'onyamalimba.png' },
			balance: 33450,
		},
		{
			id: 1004,
			name: 'Donette Foller',
			country: { name: 'Germany', code: 'de' },
			company: 'Printing Dimensions',
			date: '2021-03-18',
			status: 'renewal',
			verified: true,
			activity: 9,
			representative: { name: 'Xuxue Feng', image: 'xuxuefeng.png' },
			balance: 92100,
		},
		{
			id: 1005,
			name: 'Simona Morasca',
			country: { name: 'France', code: 'fr' },
			company: 'Chapman, Ross E Esq',
			date: '2017-08-30',
			status: 'qualified',
			verified: true,
			activity: 61,
			representative: { name: 'Stephen Shaw', image: 'stephenshaw.png' },
			balance: 13200,
		},
		{
			id: 1006,
			name: 'Mitsue Tollner',
			country: { name: 'India', code: 'in' },
			company: 'Morlong Associates',
			date: '2016-12-25',
			status: 'unqualified',
			verified: false,
			activity: 28,
			representative: { name: 'Ioni Bowcher', image: 'ionibowcher.png' },
			balance: 38400,
		},
		{
			id: 1007,
			name: 'Leota Dilliard',
			country: { name: 'Italy', code: 'it' },
			company: 'Commercial Press',
			date: '2022-01-10',
			status: 'negotiation',
			verified: true,
			activity: 39,
			representative: { name: 'Xuxue Feng', image: 'xuxuefeng.png' },
			balance: 74500,
		},
		{
			id: 1008,
			name: 'Sage Wieser',
			country: { name: 'Japan', code: 'jp' },
			company: 'Truhlar And Truhlar Attys',
			date: '2023-06-03',
			status: 'renewal',
			verified: false,
			activity: 7,
			representative: { name: 'Asiya Javayant', image: 'asiyajavayant.png' },
			balance: 6200,
		},
		{
			id: 1009,
			name: 'Kris Marrier',
			country: { name: 'South Korea', code: 'kr' },
			company: 'King, Christopher A Esq',
			date: '2022-09-14',
			status: 'qualified',
			verified: true,
			activity: 66,
			representative: { name: 'Amy Elsner', image: 'amyelsner.png' },
			balance: 81230,
		},

		{
			id: 1000,
			name: 'James Butt',
			country: { name: 'Algeria', code: 'dz' },
			company: 'Benton, John B Jr',
			date: '2015-09-13',
			status: 'unqualified',
			verified: true,
			activity: 17,
			representative: { name: 'Ioni Bowcher', image: 'ionibowcher.png' },
			balance: 70663,
		},
		{
			id: 1001,
			name: 'Josephine Darakjy',
			country: { name: 'Argentina', code: 'ar' },
			company: 'Chanay, Jeffrey A Esq',
			date: '2019-01-05',
			status: 'qualified',
			verified: false,
			activity: 42,
			representative: { name: 'Amy Elsner', image: 'amyelsner.png' },
			balance: 45890,
		},
		{
			id: 1002,
			name: 'Art Venere',
			country: { name: 'Brazil', code: 'br' },
			company: 'Chemel, James L Cpa',
			date: '2020-07-21',
			status: 'new',
			verified: true,
			activity: 12,
			representative: { name: 'Asiya Javayant', image: 'asiyajavayant.png' },
			balance: 11890,
		},
		{
			id: 1003,
			name: 'Lenna Paprocki',
			country: { name: 'Canada', code: 'ca' },
			company: 'Feltz Printing Service',
			date: '2018-11-02',
			status: 'negotiation',
			verified: false,
			activity: 54,
			representative: { name: 'Onyama Limba', image: 'onyamalimba.png' },
			balance: 33450,
		},
		{
			id: 1004,
			name: 'Donette Foller',
			country: { name: 'Germany', code: 'de' },
			company: 'Printing Dimensions',
			date: '2021-03-18',
			status: 'renewal',
			verified: true,
			activity: 9,
			representative: { name: 'Xuxue Feng', image: 'xuxuefeng.png' },
			balance: 92100,
		},
		{
			id: 1005,
			name: 'Simona Morasca',
			country: { name: 'France', code: 'fr' },
			company: 'Chapman, Ross E Esq',
			date: '2017-08-30',
			status: 'qualified',
			verified: true,
			activity: 61,
			representative: { name: 'Stephen Shaw', image: 'stephenshaw.png' },
			balance: 13200,
		},
		{
			id: 1006,
			name: 'Mitsue Tollner',
			country: { name: 'India', code: 'in' },
			company: 'Morlong Associates',
			date: '2016-12-25',
			status: 'unqualified',
			verified: false,
			activity: 28,
			representative: { name: 'Ioni Bowcher', image: 'ionibowcher.png' },
			balance: 38400,
		},
		{
			id: 1007,
			name: 'Leota Dilliard',
			country: { name: 'Italy', code: 'it' },
			company: 'Commercial Press',
			date: '2022-01-10',
			status: 'negotiation',
			verified: true,
			activity: 39,
			representative: { name: 'Xuxue Feng', image: 'xuxuefeng.png' },
			balance: 74500,
		},
		{
			id: 1008,
			name: 'Sage Wieser',
			country: { name: 'Japan', code: 'jp' },
			company: 'Truhlar And Truhlar Attys',
			date: '2023-06-03',
			status: 'renewal',
			verified: false,
			activity: 7,
			representative: { name: 'Asiya Javayant', image: 'asiyajavayant.png' },
			balance: 6200,
		},
		{
			id: 1009,
			name: 'Kris Marrier',
			country: { name: 'South Korea', code: 'kr' },
			company: 'King, Christopher A Esq',
			date: '2022-09-14',
			status: 'qualified',
			verified: true,
			activity: 66,
			representative: { name: 'Amy Elsner', image: 'amyelsner.png' },
			balance: 81230,
		},
	]

	loading.value = false
})

const getSeverity = (status) => {
	switch (status) {
		case 'unqualified':
			return 'danger'

		case 'qualified':
			return 'success'

		case 'new':
			return 'info'

		case 'negotiation':
			return 'warn'

		case 'renewal':
			return null
	}
}
</script>
