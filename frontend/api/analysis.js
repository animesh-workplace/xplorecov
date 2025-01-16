export function useUserAnalysis() {
	const config = useRuntimeConfig()
	const BASEURL = `${config.public.API_BASE_URL}`

	const uploadAnalysis = async (formData) => {
		try {
			const csrfToken = useCookie('csrftoken')

			const { data, error } = await useFetch(`${BASEURL}/job/run-main-workflow/`, {
				method: 'POST',
				body: formData,
				headers: { 'X-CSRFToken': csrfToken.value },
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

	const getAnalysis = async () => {
		try {
			const csrfToken = useCookie('csrftoken')

			const { data, error } = await useFetch(
				`${BASEURL}/job/get-submitted-workflow/?user_id=${useCookie('session').value}`,
				{
					method: 'GET',
					headers: { 'X-CSRFToken': csrfToken.value },
				},
			)

			if (error.value) {
				throw new Error(error.value || 'An error occurred')
			}

			return data.value
		} catch (err) {
			console.error(err)
			throw err
		}
	}

	const getSpecificAnalysis = async (analysis_id) => {
		try {
			const csrfToken = useCookie('csrftoken')

			const { data, error } = await useFetch(
				`${BASEURL}/job/get-specific-workflow/?user_id=${useCookie('session').value}&analysis_id=${analysis_id}`,
				{
					method: 'GET',
					headers: { 'X-CSRFToken': csrfToken.value },
				},
			)

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
		getAnalysis,
		uploadAnalysis,
		getSpecificAnalysis,
	}
}
