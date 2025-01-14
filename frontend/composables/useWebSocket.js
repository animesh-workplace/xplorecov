const activeConnections = new Map()

export const useWebSocket = (url) => {
	const error = ref(null)
	const socket = ref(null)
	const messages = ref([])
	const isConnected = ref(false)
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
		// Check if there's an existing active connection for this URL
		if (activeConnections.has(url)) {
			const existingSocket = activeConnections.get(url)
			// Check if the existing connection is still open
			if (
				existingSocket.readyState === WebSocket.OPEN ||
				existingSocket.readyState === WebSocket.CONNECTING
			) {
				console.log('Reusing existing WebSocket connection')
				socket.value = existingSocket
				isConnected.value = existingSocket.readyState === WebSocket.OPEN
				return
			} else {
				// Clean up the dead connection
				activeConnections.delete(url)
			}
		}

		try {
			socket.value = new WebSocket(url)
			activeConnections.set(url, socket.value)

			socket.value.onopen = () => {
				isConnected.value = true
				error.value = null
				console.log('WebSocket connected')
			}

			socket.value.onmessage = (event) => {
				const message = JSON.parse(event.data)
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
				messages.value.push(message)
			}

			socket.value.onclose = () => {
				isConnected.value = false
				console.log('WebSocket disconnected')
				if (!proper_disconnection) {
					setTimeout(connect, 5000)
				}
			}

			socket.value.onerror = (event) => {
				error.value = 'WebSocket error occurred'
				console.error('WebSocket error:', event)
			}
		} catch (err) {
			error.value = 'Failed to connect to WebSocket'
			console.error('Connection error:', err)
		}
	}

	const disconnect = () => {
		if (socket.value) {
			proper_disconnection.value = true
			socket.value.close()
			socket.value = null
		}
	}

	onMounted(() => {
		connect()
	})

	onUnmounted(() => {
		disconnect()
	})

	return {
		isConnected,
		messages,
		error,
		connect,
		disconnect,
		analysis_steps,
	}
}
