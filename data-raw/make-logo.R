# data-raw/make-logo.R
#
# psAve hex sticker
# 素材: data-raw/braid.png（3本の糸が撚り合わさって1本の紐になるイラスト）
# 出力: man/figures/logo.png（README・pkgdown 用）
#
# パッケージルートで実行すること: Rscript data-raw/make-logo.R

required_packages <- c("ggplot2", "hexSticker", "magick")

not_installed <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(not_installed) > 0L) {
  install.packages(not_installed)
}

library(hexSticker)
library(ggplot2)
library(magick)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

border_color <- "#F2EAD8"   # 生成り（イラストの糸色に合わせる）
text_color <- "#F2EAD8"
url_color <- "#DDE5DA"

# ------------------------------------------------------------------
# 1. 素材画像から編み紐を中心に正方形を切り出す
# ------------------------------------------------------------------

src <- image_read("data-raw/braid.png")
sq_path <- file.path(tempdir(), "braid_square.png")
image_write(image_crop(src, "1024x1024+230+0"), sq_path)

# ------------------------------------------------------------------
# 2. hexSticker 自身に「六角形だけ」を描かせるヘルパー
#    （本体レンダリングと形状が正確に一致するマスク・枠線を得る）
# ------------------------------------------------------------------

empty_plot <- ggplot() + theme_void() + theme_transparent()

render_hex_only <- function(fill, colr, h_size) {
  path <- tempfile(fileext = ".png")
  sticker(
    empty_plot, package = "",
    s_x = 1, s_y = 1, s_width = 0.1, s_height = 0.1,
    h_fill = fill, h_color = colr, h_size = h_size,
    white_around_sticker = FALSE,
    filename = path, dpi = 300
  )
  image_read(path)
}

# ------------------------------------------------------------------
# 3. 本体: 画像を六角形いっぱいに敷き，マスクで切り抜き，枠線を重ねる
# ------------------------------------------------------------------

tmp <- tempfile(fileext = ".png")
sticker(
  sq_path,
  package = "psAve",
  p_x = 1.32, p_y = 1.38, p_size = 20,
  p_color = text_color, p_family = "Aller_Rg", p_fontface = "bold",
  s_x = 1.00, s_y = 1.00, s_width = 1.20, s_height = 1.20,
  h_fill = "#31402F", h_color = NA,
  url = "kabajiro.github.io/psAve",
  u_x = 0.52, u_y = 0.40, u_size = 4.2,
  u_color = url_color, u_family = "Aller_Rg", u_angle = -30,
  white_around_sticker = FALSE,
  filename = tmp, dpi = 300
)
img <- image_read(tmp)

# マスク: 塗りつぶし六角形を白背景に落として反転（内側=白）
mask <- render_hex_only("black", NA, 1.2)
mask_gray <- image_negate(image_flatten(image_background(mask, "white")))
mask_gray <- image_convert(mask_gray, type = "grayscale", matte = FALSE)
img <- image_composite(img, mask_gray, operator = "copyopacity")

# 枠線: 縁だけの六角形を上に重ねる
ring <- render_hex_only(NA, border_color, 1.8)
img <- image_composite(img, ring, operator = "over")

image_write(img, "man/figures/logo.png")
cat("wrote man/figures/logo.png\n")
