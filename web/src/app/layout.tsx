import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Recessive Variant Explorer: non-additive effects in cattle",
  description:
    "An independent scrollytelling walkthrough of the Reynolds et al. (2021) Nature Genetics non-additive cattle GWAS, with reproduction, verification and a mating-risk simulator.",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en-NZ">
      <head>
        <script
          dangerouslySetInnerHTML={{
            __html:
              "try{var t=localStorage.getItem('theme');if(t){document.documentElement.setAttribute('data-theme',t);}}catch(e){}",
          }}
        />
      </head>
      <body>{children}</body>
    </html>
  );
}
