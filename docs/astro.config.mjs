// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

export default defineConfig({
  site: 'https://psy-fer.github.io',
  base: '/rustar-aligner',
  integrations: [
    starlight({
      title: 'rustar-aligner',
      logo: {
        src: './src/assets/rustar-icon.svg',
      },
      favicon: '/favicon.svg',
      components: {
        Header: './src/components/Header.astro',
      },
      customCss: ['./src/styles/custom.css'],
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/Psy-Fer/rustar-aligner',
        },
      ],
      editLink: {
        baseUrl: 'https://github.com/Psy-Fer/rustar-aligner/edit/main/docs/',
      },
      sidebar: [
        {
          label: 'Getting started',
          items: [
            { label: 'Introduction', slug: 'getting-started/introduction' },
            { label: 'Installation', slug: 'getting-started/installation' },
            { label: 'Quick start', slug: 'getting-started/quick-start' },
          ],
        },
        {
          label: 'Guides',
          items: [
            { label: 'Generating a genome index', slug: 'guides/genome-index' },
            { label: 'Single-end alignment', slug: 'guides/single-end' },
            { label: 'Paired-end alignment', slug: 'guides/paired-end' },
            { label: 'Two-pass mode', slug: 'guides/two-pass' },
            { label: 'Chimeric detection', slug: 'guides/chimeric' },
            { label: 'Gene quantification', slug: 'guides/quantification' },
            { label: 'Migrating from STAR', slug: 'guides/migrating-from-star' },
          ],
        },
        {
          label: 'Reference',
          items: [
            { label: 'CLI parameters', slug: 'reference/cli-parameters' },
            { label: 'Output files', slug: 'reference/output-files' },
            { label: 'STAR compatibility', slug: 'reference/star-compatibility' },
            { label: 'Performance', slug: 'reference/performance' },
          ],
        },
        {
          label: 'About',
          items: [
            { label: 'Contributing', slug: 'about/contributing' },
            { label: 'Changelog', slug: 'about/changelog' },
            { label: 'License', slug: 'about/license' },
          ],
        },
      ],
    }),
  ],
});
