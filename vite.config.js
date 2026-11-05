import { defineConfig } from 'vite';
import { viteSingleFile } from 'vite-plugin-singlefile';

export default defineConfig({
  // `npm run build` inlines everything (JS, CSS, samples) into dist/index.html,
  // so the app opens straight from disk with a double-click, no server needed.
  plugins: [viteSingleFile()],
  server: { port: 5173 },
  build: { chunkSizeWarningLimit: 1000 }, // three.js alone is ~600 kB; that's expected
});
