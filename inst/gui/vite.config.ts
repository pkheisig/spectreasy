import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

const configuredPort = Number.parseInt(process.env.VITE_DEV_PORT || '5174', 10)

// https://vite.dev/config/
export default defineConfig({
  base: process.env.VITE_BASE_PATH || '/',
  plugins: [
    react(),
    tailwindcss(),
  ],
  server: {
    host: '127.0.0.1',
    port: Number.isFinite(configuredPort) ? configuredPort : 5174,
    strictPort: true,
    watch: {
      usePolling: true,
    },
  },
})
