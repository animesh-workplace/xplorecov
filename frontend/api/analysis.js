export function useUserAnalysis() {
	const config = useRuntimeConfig()
	const BASEURL = `${config.public.API_BASE_URL}`

	const uploadAnalysis = async (formData) => {
		try {
			const { data, error } = await useFetch(`${BASEURL}/job/run-main-workflow/`, {
				method: 'POST',
				body: formData,
			})

			if (error.value) {
				throw new Error(error.value || 'An error occurred')
			}

			return data.value
		} catch (err) {
			console.error(err)
			throw err
		}
	}

	return {
		uploadAnalysis,
	}
}
