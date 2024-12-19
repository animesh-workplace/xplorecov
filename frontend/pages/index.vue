<template>
  <div>
    <Button
      @click="refreshData"
      rounded
      size="large"
      label="Large button"
      icon="pi pi-check"
    />

    <div v-if="isLoading" class="loader">Loading Graph...</div>

    <VChart
      ref="chart"
      class="test"
      :option="option"
      :class="{ hidden: isLoading }"
    />

    <ReadTable />
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, nextTick } from "vue";
const chart = templateRef("chart");

function random() {
  return Math.round(300 + Math.random() * 700) / 10;
}

function getData(): ECOption {
  return {
    animation: true,
    tooltip: {
      className: "echarts-tooltip",
    },
    toolbox: {
      show: true,
      feature: {
        dataZoom: {},
        saveAsImage: {},
      },
    },
    dataset: {
      dimensions: ["Product", "2015", "2016", "2017"],
      source: [
        {
          Product: "Matcha Latte",
          2015: random(),
          2016: random(),
          2017: random(),
        },
        {
          Product: "Milk Tea",
          2015: random(),
          2016: random(),
          2017: random(),
        },
        {
          Product: "Cheese Cocoa",
          2015: random(),
          2016: random(),
          2017: random(),
        },
        {
          Product: "Walnut Brownie",
          2015: random(),
          2016: random(),
          2017: random(),
        },
      ],
    },
    xAxis: { type: "category" },
    yAxis: {},
    itemStyle: { borderRadius: [5, 5, 0, 0] },
    series: [{ type: "bar" }, { type: "bar" }, { type: "bar" }],
  };
}

const isLoading = ref(true);

const option = shallowRef(getData());
function refreshData() {
  option.value = getData();
}

onMounted(() => {
  nextTick(() => {
    isLoading.value = false;
  });
});
</script>

<style>
.test {
  margin: 5px;
  width: 700px;
  height: 500px;
}
</style>
