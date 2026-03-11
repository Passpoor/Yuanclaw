# WWW文件夹 - Shiny静态资源

## 文件夹说明
此文件夹是Shiny应用的标准静态文件目录，用于存放可被Web浏览器直接访问的静态资源。

## 功能
- Shiny自动将www文件夹中的文件映射到应用根路径
- 可以通过相对路径直接访问：`logo.png`
- 支持所有Web静态资源：图片、CSS、JS等

## 当前文件
- `logo.png` - YuanSeq 应用 Logo

## 添加其他静态资源
如需添加其他静态文件（如CSS、JS、图片等），可以直接放入此文件夹。

## 路径规则
- HTML中引用：`<img src="logo.png">`
- Shiny中引用：`tags$img(src="logo.png")`
- 自动映射：`www/logo.png` → `http://localhost:port/logo.png`

## 注意事项
- 此文件夹中的所有文件都是公开可访问的
- 不要放置敏感信息在此文件夹中
- 文件名建议使用小写字母和连字符