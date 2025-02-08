export function useTyping(text, speed = 100) {
	const displayedText = ref('')

	watchEffect(async () => {
		displayedText.value = '' // Reset text
		for (let i = 0; i < text.length; i++) {
			await new Promise((resolve) => setTimeout(resolve, speed))
			displayedText.value += text[i]
		}
	})

	return { displayedText }
}
