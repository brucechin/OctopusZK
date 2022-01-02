
for TOTAL_CORES in 16; do
  for SIZE in `seq 15 15`; do

    export APP=fmsm-g1
    export MEMORY=64G
    export MULTIPLIER=1

    export CORES=1
    export NUM_EXECUTORS=$((TOTAL_CORES / CORES))
    export NUM_PARTITIONS=$((TOTAL_CORES * MULTIPLIER))

    sudo /opt/spark/bin/spark-submit \
      --conf spark.driver.memory=$MEMORY \
      --conf spark.driver.maxResultSize=$MEMORY \
      --conf spark.executor.cores=$CORES \
      --total-executor-cores $TOTAL_CORES \
      --conf spark.executor.memory=$MEMORY \
      --conf spark.memory.fraction=0.95 \
      --conf spark.memory.storageFraction=0.3 \
      --conf spark.kryoserializer.buffer.max=1g \
      --conf spark.rdd.compress=true \
      --conf spark.rpc.message.maxSize=2000 \
      --conf spark.executor.heartbeatInterval=30s \
      --conf spark.network.timeout=300s\
      --conf spark.speculation=true \
      --conf spark.speculation.interval=5000ms \
      --conf spark.speculation.multiplier=2 \
      --conf spark.local.dir=/mnt/spark \
      --conf spark.logConf=true \
      --conf spark.eventLog.dir=/tmp/spark-events \
      --conf spark.eventLog.enabled=false \
      --class "profiler.Profiler" \
      ./target/dizk.jar $NUM_EXECUTORS $CORES $MEMORY $APP $SIZE $NUM_PARTITIONS
  done
done
