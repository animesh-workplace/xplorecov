import { map } from 'lodash'

export const useWebSocket = (url) => {
	const socket = ref(null)
	const tools_version = ref({ tools: {} })
	const proper_disconnection = ref(false)
	const analysis_steps = ref({
		// Status - Pending, Loading, Completed
		step1: { index: 1, name: 'Analysis Queued', status: 'pending' },
		step2: { index: 2, name: 'Runnning QC Checks', status: 'pending' },
		step3: { index: 3, name: 'Updating Tools', status: 'pending' },
		step4: { index: 4, name: 'Runnning Nextclade Analysis', status: 'pending' },
		step5: { index: 5, name: 'Running Pangolin Analysis', status: 'pending' },
		step6: { index: 6, name: 'Summarizing Results', status: 'pending' },
	})

	const connect = () => {
		try {
			socket.value = new WebSocket(url)

			socket.value.onopen = () => {
				console.log('WebSocket connected')
			}

			socket.value.onmessage = (event) => {
				const message = JSON.parse(event.data)
				if ('tools_used' in message) {
					tools_version.value = message.tools_used
				}
				if ('current_status' in message) {
					map(message.current_status, (d) => update_analysis_steps({ message: d }))
				} else if (typeof message.message == 'object') {
					update_analysis_steps(message)
				}
			}

			socket.value.onclose = () => {
				console.log('WebSocket disconnected')
				if (!proper_disconnection.value) {
					setTimeout(connect, 5000)
				}
			}

			socket.value.onerror = (event) => {
				console.error('WebSocket error:', event)
			}
		} catch (err) {
			console.error('Connection error:', err)
		}
	}

	const update_analysis_steps = (message) => {
		let update_step = null
		if (typeof message.message == 'object') {
			if (message.message.type == 'nextclade-rule') {
				update_step = 'step4'
			} else if (message.message.type == 'pangolin-usher') {
				update_step = 'step5'
			} else if (message.message.type == 'combine') {
				update_step = 'step6'
			} else if (message.message.type == 'WORKFLOW') {
				if (message.message.status == 'Started Workflow') {
					analysis_steps.value.step1.status = 'completed'
					analysis_steps.value.step2.status = 'completed'
					analysis_steps.value.step3.status = 'completed'
				}
			}

			if (message.message.status == 'start') {
				analysis_steps.value[update_step].status = 'loading'
			} else if (message.message.status == 'end') {
				analysis_steps.value[update_step].status = 'completed'
			}
		}
	}

	const disconnect = () => {
		if (socket.value) {
			socket.value.close()
			socket.value = null
		}
	}

	onMounted(() => {
		connect()
	})

	onUnmounted(() => {
		proper_disconnection.value = true
		disconnect()
	})

	return {
		analysis_steps,
		tools_version,
	}
}
