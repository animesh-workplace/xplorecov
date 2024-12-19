import Aura from "@primevue/themes/aura";

export default defineNuxtConfig({
  compatibilityDate: "2024-11-01",
  devtools: { enabled: true },
  modules: [
    "@primevue/nuxt-module",
    "nuxt-echarts",
    "@pinia/nuxt",
    "@nuxtjs/tailwindcss",
    "@nuxt/image",
    "@nuxt/icon",
    "@vueuse/nuxt",
    "nuxt-umami",
    "dayjs-nuxt",
  ],
  primevue: {
    options: {
      ripple: true,
      theme: { preset: Aura },
    },
  },
  umami: { enabled: false },
  echarts: {
    ssr: true,
    charts: ["BarChart"],
    components: [
      "GridComponent",
      "DatasetComponent",
      "TooltipComponent",
      "ToolboxComponent",
    ],
  },
});
