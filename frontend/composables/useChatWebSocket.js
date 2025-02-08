export const useChatWebSocket = (url) => {
	const message = ref({ message: '' })
	const socket = ref(null)
	const proper_disconnection = ref(false)

	const connect = () => {
		try {
			socket.value = new WebSocket(url)

			socket.value.onopen = () => {
				console.log('WebSocket connected')
			}

			socket.value.onmessage = (event) => {
				message.value = JSON.parse(event.data)
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

	const Chatdisconnect = (value_for_disconnection) => {
		if (socket.value) {
			socket.value.close()
			socket.value = null
			if (value_for_disconnection) {
				proper_disconnection.value = value_for_disconnection
			}
		}
	}

	onMounted(() => {
		connect()
	})

	onUnmounted(() => {
		proper_disconnection.value = true
		Chatdisconnect()
	})

	return {
		Chatdisconnect,
		message,
	}
}
