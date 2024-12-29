import { defineStore } from 'pinia'
import { v4 as uuidv4 } from 'uuid'

export const useSessionStore = defineStore('session', {
	state: () => ({ session: undefined }),
	actions: {
		setSessionInfo() {
			this.session = uuidv4()
			const session_cookie = useCookie('session')
			session_cookie.value = this.session
		},
		getSessionInfo() {
			const session_cookie = useCookie('session')
			if (session_cookie.value) {
				this.session = session_cookie.value
			} else {
				this.setSessionInfo()
			}
		},
	},
})
