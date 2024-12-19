export function usePondFileError() {
	/**
	 * Handles setting an error or warning state for a FilePond file item.
	 * @param {HTMLElement} reference - The FilePond instance reference. This is a required parameter.
	 * @param {string} type - The type of state (options are - 'processing-error', 'processing-warn', 'processing-complete'). This is a required parameter.
	 * @param {string} main_text - The main message text.
	 * @param {string} sub_text - The sub-message text.
	 *
	 * @throws {Error} Throws an error if the reference or type is missing or invalid.
	 */

	const handlePondFileError = (reference, type, mainText = '', subText = '') => {
		// Check if required parameters are provided
		if (!reference || !(reference?.$el instanceof HTMLElement)) {
			throw new Error('The "reference" parameter must be a valid HTMLElement.')
		}

		const validTypes = ['processing-error', 'processing-warn', 'processing-complete']
		if (!type || !validTypes.includes(type)) {
			throw new Error(
				'The "type" parameter must be one of: processing-error, processing-warn, processing-complete.',
			)
		}

		const fileItem = reference.$el.querySelector('.filepond--item')
		const status = reference.$el.querySelector('.filepond--file-status')
		const statusMain = reference.$el.querySelector('.filepond--file-status-main')
		const statusSub = reference.$el.querySelector('.filepond--file-status-sub')

		// Set the file item state
		if (fileItem) fileItem.setAttribute('data-filepond-item-state', type)

		// Update file status visibility
		if (status) {
			status.style.opacity = 1
			status.style.visibility = 'visible'
		}

		// Update sub status opacity
		if (statusSub) statusSub.style.opacity = 1

		// Update status text
		if (statusMain) statusMain.textContent = mainText
		if (statusSub) statusSub.textContent = subText
	}

	return { handlePondFileError }
}
