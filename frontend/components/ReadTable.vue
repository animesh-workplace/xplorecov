<template>
  <div class="mx-20">
    <div v-if="isLoading" class="loader">Loading...</div>

    <div>
      <file-pond
        :class="{ hidden: isLoading }"
        ref="pond"
        credits="false"
        allow-multiple="true"
        v-model:files="myFiles"
        @init="handleFilePondInit"
        label-idle="Drop files here..."
        accepted-file-types="image/jpeg, image/png"
        labelIdle="
            <span class='is-family-primary has-text-weight-semibold has-text-grey-dark is-clickable'>
                Drag & Drop your Metadata or
            </span>
            <span
                class='is-family-primary has-text-weight-semibold has-text-grey-dark is-clickable'
                style='text-decoration: underline;'
            >
                Browse
            </span>
		"
      />
    </div>
  </div>
</template>

<script setup>
// Import FilePond
import vueFilePond from "vue-filepond";

// Import plugins
import FilePondPluginFileValidateType from "filepond-plugin-file-validate-type/dist/filepond-plugin-file-validate-type.esm.js";
import FilePondPluginImagePreview from "filepond-plugin-image-preview/dist/filepond-plugin-image-preview.esm.js";

// Import styles
import "filepond/dist/filepond.min.css";
import "filepond-plugin-image-preview/dist/filepond-plugin-image-preview.min.css";

// Create FilePond component
const FilePond = vueFilePond(
  FilePondPluginFileValidateType,
  FilePondPluginImagePreview
);

// Reactive state for files
const myFiles = ref([]);
const isLoading = ref(true);

// Methods
const handleFilePondInit = () => {
  console.log("FilePond has initialized");

  isLoading.value = false;
  // example of instance method call on pond reference
  pondRef.value?.getFiles();
};

// Reference to FilePond component
const pondRef = ref(null);
</script>
