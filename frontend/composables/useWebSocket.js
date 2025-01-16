import { map, round } from 'lodash-es'

export const useWebSocket = (url) => {
	const socket = ref(null)
	const dayjs = useDayjs()
	const tools_version = ref({ tools: {} })
	const proper_disconnection = ref(false)
	const analysis_steps = ref({
		// Status - Pending, Loading, Completed
		step1: { index: 1, name: 'Quality Control Checks', status: 'completed', duration: 0 },
		step2: { index: 2, name: 'Tool Updates', status: 'completed', duration: 0 },
		step3: { index: 3, name: 'Queuing Analysis', status: 'pending', duration: 0 },
		step4: { index: 4, name: 'Nextclade Analysis Execution', status: 'pending', duration: 0 },
		step5: { index: 5, name: 'Pangolin Analysis Execution', status: 'pending', duration: 0 },
		step6: { index: 6, name: 'Results Summarization', status: 'pending', duration: 0 },
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
		if (typeof message.message == 'object') {
			analysis_steps.value[message.message.step_id].status = 'loading'
			if (message.message.status == 'start') {
				analysis_steps.value[message.message.step_id].status = 'loading'
				analysis_steps.value[message.message.step_id].duration = dayjs(message.message.timestamp)
			} else if (message.message.status == 'end') {
				analysis_steps.value[message.message.step_id].status = 'completed'
				const durationInSeconds = round(
					dayjs(message.message.timestamp).diff(
						analysis_steps.value[message.message.step_id].duration,
						'seconds',
						true,
					),
					2,
				)
				analysis_steps.value[message.message.step_id].duration =
					durationInSeconds < 60
						? durationInSeconds
						: round(
								dayjs(message.message.timestamp).diff(
									analysis_steps.value[message.message.step_id].duration,
									'minutes',
									true,
								),
								2,
							)
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
