import type { NextConfig } from "next";

// Static export: the site is plain files (no server, no API keys, no runtime LLM calls).
const nextConfig: NextConfig = {
  output: "export",
  images: { unoptimized: true },
  reactStrictMode: true,
};

export default nextConfig;
