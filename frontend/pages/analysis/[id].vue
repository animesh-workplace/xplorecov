<template>
	<section>
		<div class="bg-surface-0 dark:bg-surface-950 px-6 py-8 md:px-12 lg:px-20">
			<ul class="list-none p-0 m-0 flex items-center font-medium mb-4">
				<li>
					<a class="text-surface-500 dark:text-surface-300 no-underline leading-normal cursor-pointer"
						>Application</a
					>
				</li>
				<li class="px-2">
					<i class="pi pi-angle-right text-surface-500 dark:text-surface-300 leading-normal" />
				</li>
				<li>
					<span class="text-surface-900 dark:text-surface-0 leading-normal">Analytics</span>
				</li>
			</ul>
			<div class="flex items-start flex-col lg:justify-between lg:flex-row">
				<div>
					<div class="font-medium text-3xl text-surface-900 dark:text-surface-0">Customers</div>
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
				<!-- <div class="mt-4 lg:mt-0">
					<Button label="Add" class="mr-2" outlined icon="pi pi-user-plus" />
					<Button label="Save" icon="pi pi-check" />
				</div> -->
			</div>
		</div>

		<div class="flex gap-4 justify-center px-6 py-8 md:px-12 lg:px-20">
			<div
				:key="item"
				v-for="item in [0, 1, 2, 3, 4, 5]"
				class="relative w-64 p-4 overflow-hidden bg-white shadow-lg rounded-xl"
			>
				<img
					alt="image"
					src="https://www.tailwind-kit.com/images/object/1.png"
					class="absolute w-40 h-40 mb-4 -right-16 -bottom-16"
				/>
				<div class="w-4/6">
					<p class="mb-2 text-lg font-medium text-gray-800">NextJS</p>
					<p class="text-xs text-gray-400">NextJs build all free</p>
				</div>
			</div>
		</div>

		<div class="px-6 py-8 md:px-12 lg:px-20">
			<div class="w-full bg-white shadow-lg rounded-2xl dark:bg-gray-700">
				<p class="p-4 font-bold text-black text-md dark:text-white">My Tasks</p>
				<ul>
					<li
						:key="index"
						v-for="(step, index) in analysis_steps"
						class="flex items-center justify-between py-3 text-gray-600 border-gray-100 dark:text-gray-200 dark:border-gray-800"
						:class="{ 'border-b-2': step.index != 6, 'border-t-2': step.index == 1 }"
					>
						<div class="flex items-center justify-start text-sm">
							<span class="mx-4"> {{ step.index }} </span>
							<span> {{ step.name }} </span>
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

		<!-- <div class="px-6 py-8 md:px-12 lg:px-20">
			<DataTable
				paginator
				:rows="5"
				dataKey="id"
				showGridlines
				striped-rows
				removableSort
				:value="customers"
				:loading="loading"
				class="!rounded-lg"
				filterDisplay="menu"
				v-model:filters="filters"
				v-model:expandedRows="expandedRows"
				:globalFilterFields="['name', 'country.name', 'representative.name', 'balance', 'status']"
			>
				<template #header>
					<div class="flex justify-between">
						<Button
							type="button"
							icon="pi pi-filter-slash"
							label="Clear"
							outlined
							@click="clearFilter()"
						/>
						<IconField>
							<InputIcon>
								<i class="pi pi-search" />
							</InputIcon>
							<InputText v-model="filters['global'].value" placeholder="Keyword Search" />
						</IconField>
					</div>
				</template>
				<Column field="name" header="Name" style="min-width: 12rem">
					<template #body="{ data }">
						{{ data.name }}
					</template>
					<template #filter="{ filterModel }">
						<InputText v-model="filterModel.value" type="text" placeholder="Search by name" />
					</template>
				</Column>

				<Column header="Country" filterField="country.name" style="min-width: 12rem">
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
					<template #filter="{ filterModel }">
						<InputText v-model="filterModel.value" type="text" placeholder="Search by country" />
					</template>
					<template #filterclear="{ filterCallback }">
						<Button
							type="button"
							icon="pi pi-times"
							@click="filterCallback()"
							severity="secondary"
						></Button>
					</template>
					<template #filterapply="{ filterCallback }">
						<Button
							type="button"
							icon="pi pi-check"
							@click="filterCallback()"
							severity="success"
						></Button>
					</template>
					<template #filterfooter>
						<div class="px-4 pt-0 pb-4 text-center">Customized Buttons</div>
					</template>
				</Column>

				<Column
					header="Agent"
					style="min-width: 14rem"
					filterField="representative"
					:showFilterMatchModes="false"
					:filterMenuStyle="{ width: '14rem' }"
				>
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
					<template #filter="{ filterModel }">
						<MultiSelect
							v-model="filterModel.value"
							:options="representatives"
							optionLabel="name"
							placeholder="Any"
						>
							<template #option="slotProps">
								<div class="flex items-center gap-2">
									<img
										:alt="slotProps.option.name"
										:src="`https://primefaces.org/cdn/primevue/images/avatar/${slotProps.option.image}`"
										style="width: 32px"
									/>
									<span>{{ slotProps.option.name }}</span>
								</div>
							</template>
						</MultiSelect>
					</template>
				</Column>

				<Column header="Date" filterField="date" dataType="date" style="min-width: 10rem">
					<template #body="{ data }">
						{{ formatDate(data.date) }}
					</template>
					<template #filter="{ filterModel }">
						<DatePicker v-model="filterModel.value" dateFormat="mm/dd/yy" placeholder="mm/dd/yyyy" />
					</template>
				</Column>

				<Column header="Balance" filterField="balance" dataType="numeric" style="min-width: 10rem">
					<template #body="{ data }">
						{{ formatCurrency(data.balance) }}
					</template>
					<template #filter="{ filterModel }">
						<InputNumber v-model="filterModel.value" mode="currency" currency="USD" locale="en-US" />
					</template>
				</Column>

				<Column
					header="Status"
					field="status"
					:filterMenuStyle="{ width: '14rem' }"
					style="min-width: 12rem"
				>
					<template #body="{ data }">
						<Tag :value="data.status" :severity="getSeverity(data.status)" />
					</template>
					<template #filter="{ filterModel }">
						<Select v-model="filterModel.value" :options="statuses" placeholder="Select One" showClear>
							<template #option="slotProps">
								<Tag :value="slotProps.option" :severity="getSeverity(slotProps.option)" />
							</template>
						</Select>
					</template>
				</Column>

				<Column field="activity" header="Activity" :showFilterMatchModes="false" style="min-width: 12rem">
					<template #body="{ data }">
						<ProgressBar :value="data.activity" :showValue="false" style="height: 6px"></ProgressBar>
					</template>
					<template #filter="{ filterModel }">
						<Slider v-model="filterModel.value" range class="m-4"></Slider>
						<div class="flex items-center justify-between px-2">
							<span>{{ filterModel.value ? filterModel.value[0] : 0 }}</span>
							<span>{{ filterModel.value ? filterModel.value[1] : 100 }}</span>
						</div>
					</template>
				</Column>

				<Column
					field="verified"
					header="Verified"
					dataType="boolean"
					bodyClass="text-center"
					style="min-width: 8rem"
				>
					<template #body="{ data }">
						<i
							class="pi"
							:class="{
								'pi-check-circle text-green-500 ': data.verified,
								'pi-times-circle text-red-500': !data.verified,
							}"
						></i>
					</template>
					<template #filter="{ filterModel }">
						<label for="verified-filter" class="font-bold"> Verified </label>
						<Checkbox
							v-model="filterModel.value"
							:indeterminate="filterModel.value === null"
							binary
							inputId="verified-filter"
						/>
					</template>
				</Column>
			</DataTable>
		</div> -->
	</section>
</template>

<script setup>
import { ref, onMounted } from 'vue'
// import { CustomerService } from '@/service/CustomerService'
import { FilterMatchMode, FilterOperator } from '@primevue/core/api'

const customers = ref()
const filters = ref()
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
const expandedRows = ref({})
const analysis_steps = ref({
	// Status - Pending, Loading, Completed
	step1: { index: 1, name: 'Analysis Queued', status: 'pending' },
	step2: { index: 2, name: 'Runnning QC Checks', status: 'pending' },
	step3: { index: 3, name: 'Updating Tools', status: 'pending' },
	step4: { index: 4, name: 'Runnning Nextclade Analysis', status: 'pending' },
	step5: { index: 5, name: 'Running Pangolin Analysis', status: 'loading' },
	step6: { index: 6, name: 'Summarizing Results', status: 'completed' },
})

onMounted(() => {
	customers.value = [
		{
			id: 1000,
			name: 'James Butt',
			country: {
				name: 'Algeria',
				code: 'dz',
			},
			company: 'Benton, John B Jr',
			date: '2018-06-15',
			status: 'unqualified',
			verified: true,
			activity: 45,
			representative: {
				name: 'Ioni Bowcher',
				image: 'ionibowcher.png',
			},
			balance: 70663,
		},
		{
			id: 1001,
			name: 'Jane Doe',
			country: {
				name: 'United States',
				code: 'us',
			},
			company: 'Wayne Enterprises',
			date: '2020-12-11',
			status: 'qualified',
			verified: false,
			activity: 78,
			representative: {
				name: 'Amy Elsner',
				image: 'amyelsner.png',
			},
			balance: 92038,
		},
		{
			id: 1002,
			name: 'Chris Evans',
			country: {
				name: 'Brazil',
				code: 'br',
			},
			company: 'Stark Industries',
			date: '2019-03-25',
			status: 'new',
			verified: true,
			activity: 62,
			representative: {
				name: 'Anna Fali',
				image: 'annafali.png',
			},
			balance: 56022,
		},
		{
			id: 1003,
			name: 'Scarlett Johansson',
			country: {
				name: 'France',
				code: 'fr',
			},
			company: 'Oscorp',
			date: '2021-08-19',
			status: 'renewal',
			verified: false,
			activity: 34,
			representative: {
				name: 'Asiya Javayant',
				image: 'asiyajavayant.png',
			},
			balance: 48200,
		},
		{
			id: 1004,
			name: 'Tom Holland',
			country: {
				name: 'India',
				code: 'in',
			},
			company: 'Daily Bugle',
			date: '2023-01-05',
			status: 'negotiation',
			verified: true,
			activity: 19,
			representative: {
				name: 'Brad Jackson',
				image: 'asiyajavayant.png',
			},
			balance: 61200,
		},
		{
			id: 1005,
			name: 'Natalie Portman',
			country: {
				name: 'Canada',
				code: 'ca',
			},
			company: 'Acme Corp',
			date: '2016-07-22',
			status: 'unqualified',
			verified: true,
			activity: 82,
			representative: {
				name: 'Amy Elsner',
				image: 'amyelsner.png',
			},
			balance: 73500,
		},
		{
			id: 1006,
			name: 'Robert Downey Jr.',
			country: {
				name: 'Germany',
				code: 'de',
			},
			company: 'Hydra',
			date: '2017-11-09',
			status: 'qualified',
			verified: false,
			activity: 47,
			representative: {
				name: 'Anna Fali',
				image: 'annafali.png',
			},
			balance: 88400,
		},
		{
			id: 1007,
			name: 'Mark Ruffalo',
			country: {
				name: 'China',
				code: 'cn',
			},
			company: 'Pym Technologies',
			date: '2022-05-14',
			status: 'new',
			verified: true,
			activity: 51,
			representative: {
				name: 'Ioni Bowcher',
				image: 'ionibowcher.png',
			},
			balance: 49200,
		},
		{
			id: 1008,
			name: 'Zendaya',
			country: {
				name: 'Australia',
				code: 'au',
			},
			company: 'SHIELD',
			date: '2020-02-03',
			status: 'renewal',
			verified: false,
			activity: 76,
			representative: {
				name: 'Brad Jackson',
				image: 'xuxuefeng.png',
			},
			balance: 30000,
		},
		{
			id: 1009,
			name: 'John Smith',
			country: {
				name: 'Argentina',
				code: 'ar',
			},
			company: 'Parker Industries',
			date: '2015-09-13',
			status: 'negotiation',
			verified: true,
			activity: 29,
			representative: {
				name: 'Asiya Javayant',
				image: 'asiyajavayant.png',
			},
			balance: 78000,
		},
	]

	// CustomerService.getCustomersMedium().then((data) => {
	// 	customers.value = getCustomers(data)
	loading.value = false
	// })
})

const initFilters = () => {
	filters.value = {
		global: { value: null, matchMode: FilterMatchMode.CONTAINS },
		name: {
			operator: FilterOperator.AND,
			constraints: [{ value: null, matchMode: FilterMatchMode.STARTS_WITH }],
		},
		'country.name': {
			operator: FilterOperator.AND,
			constraints: [{ value: null, matchMode: FilterMatchMode.STARTS_WITH }],
		},
		representative: { value: null, matchMode: FilterMatchMode.IN },
		date: { operator: FilterOperator.AND, constraints: [{ value: null, matchMode: FilterMatchMode.DATE_IS }] },
		balance: {
			operator: FilterOperator.AND,
			constraints: [{ value: null, matchMode: FilterMatchMode.EQUALS }],
		},
		status: { operator: FilterOperator.OR, constraints: [{ value: null, matchMode: FilterMatchMode.EQUALS }] },
		activity: { value: [0, 100], matchMode: FilterMatchMode.BETWEEN },
		verified: { value: null, matchMode: FilterMatchMode.EQUALS },
	}
}

initFilters()

const formatDate = (value) => {
	return value
	return value.toLocaleDateString('en-US', {
		day: '2-digit',
		month: '2-digit',
		year: 'numeric',
	})
}
const formatCurrency = (value) => {
	return value.toLocaleString('en-US', { style: 'currency', currency: 'USD' })
}
const clearFilter = () => {
	initFilters()
}
const getCustomers = (data) => {
	return [...(data || [])].map((d) => {
		d.date = new Date(d.date)

		return d
	})
}
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

<style scoped></style>
