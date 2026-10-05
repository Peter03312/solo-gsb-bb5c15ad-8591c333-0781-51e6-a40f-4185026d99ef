# 外层目录构建入口：与 liftover-web/Dockerfile 构建同一个应用镜像。
#
# 当 Compose 项目目录为外层目录时（例如在外层目录执行：
#   docker compose --project-directory . -f liftover-web/docker-compose.yml up -d --build
# ），构建上下文是外层目录，因此 COPY 源路径带 liftover-web/ 前缀。
# 进入 liftover-web/ 后直接运行 docker compose 时，则由 liftover-web/Dockerfile
# 以该目录为上下文完成同样的构建。
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# 先装依赖，利用层缓存
COPY liftover-web/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# 拷贝应用代码与静态资源；缺少任一文件都会让构建直接失败，
# 不会产出缺少资源却显示健康的镜像
COPY liftover-web/app ./app
COPY liftover-web/static ./static

# 非 root 用户运行
RUN useradd --create-home --uid 10001 appuser \
    && chown -R appuser:appuser /app
USER appuser

EXPOSE 8080

# 容器内固定监听 8080；宿主机端口映射通过 compose / docker -p 配置
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8080/healthz', timeout=4).status==200 else 1)"

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]
