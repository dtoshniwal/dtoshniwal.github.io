# Deepesh Toshniwal Website

Personal website and blog built with [Astro](https://astro.build), based on Astro Cactus and customized for research, teaching, travel, photos, and long-form writing.

<!-- ## Overview

- Framework: Astro v6
- Styling: Tailwind CSS v4 + custom CSS
- Content: Markdown/MDX collections for posts, notes, and tags
- Search: Pagefind (generated after build)
- SEO: RSS, sitemap, robots.txt, Open Graph image generation via Satori

## Project Structure

```text
src/
  content/
    post/        # Blog posts
    note/        # Notes
    tag/         # Optional per-tag page metadata
  pages/
    about.md
    research.md
    teaching.md
    travel.md
    photos.md
    posts/
    notes/
    tags/
  components/
  layouts/
  styles/
  site.config.ts # Site metadata, menus, locale/date config
```

## Local Development

### 1. Install dependencies

```bash
pnpm install
```

### 2. Start the dev server

```bash
pnpm dev
```

### 3. Build for production

```bash
pnpm build
```

### 4. Generate local search index (Pagefind)

```bash
pnpm postbuild
```

### 5. Preview production build

```bash
pnpm preview
```

## Available Scripts

| Command | Description |
| --- | --- |
| `pnpm dev` | Run Astro in development mode |
| `pnpm start` | Alias for `pnpm dev` |
| `pnpm build` | Build static site to `dist/` |
| `pnpm postbuild` | Build Pagefind search index from `dist/` |
| `pnpm preview` | Preview built site locally |
| `pnpm check` | Type/content checks + Biome check |
| `pnpm lint` | Run Biome with `--write` |
| `pnpm format` | Format with Prettier |

## Content Authoring

This site uses Astro Content Collections defined in `src/content.config.ts`.

### Blog posts

Create Markdown/MDX files in `src/content/post/`.

Post frontmatter:

```yaml
title: "Post title"
description: "Short summary"
publishDate: 2026-01-01
updatedDate: 2026-01-15 # optional
tags: [research, fem, iga] # optional
draft: false # optional, defaults to false
pinned: false # optional, defaults to false
ogImage: "/images/custom-og.png" # optional
coverImage: # optional
  src: "./cover.png"
  alt: "Cover description"
```

### Notes

Create Markdown/MDX files in `src/content/note/`.

Note frontmatter:

```yaml
title: "Note title"
description: "Optional summary"
publishDate: 2026-01-01T00:00:00Z
```

### Tags

Optional tag metadata lives in `src/content/tag/`.
The filename is the tag slug (for example `src/content/tag/iga.md` corresponds to `/tags/iga`).

## Core Configuration

- `src/site.config.ts`
  - Site URL, title, author, description
  - Header/footer menu links
  - Date locale and formatting options
  - Expressive Code theme configuration
- `astro.config.ts`
  - Integrations and build settings
- `src/styles/global.css`
  - Global visual styles and theme variables

## Search

Search is powered by Pagefind and built from static output.

- Build first: `pnpm build`
- Then index: `pnpm postbuild`

If you skip `postbuild`, search UI loads but returns no indexed results.

## Deployment

This is a static Astro site (`output: "static"`) and deploys to any static host.

Typical flow:

1. `pnpm build`
2. `pnpm postbuild`
3. Deploy `dist/` -->

## Credits

- Original template inspiration: [Astro Cactus](https://github.com/chrismwilliams/astro-theme-cactus)

## License

MIT
